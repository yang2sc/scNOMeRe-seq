# tf enrichment

#homer size given
files<-list.files("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/homer.TF",recursive = T,pattern = "Results.txt",full.names = T)
basename(dirname(files))
foo<-lapply(files, function(f)fread(f,select = c(1,2,3)) %>% .[,anno:=basename(dirname(f)) %>% 
                                                                 sub("20190810.accrna.correlations.","",.)%>%
                                                                 sub(".high.tsv_motif_unmask_sizeGiven","",.)])%>% rbindlist()

colnames(foo)<-c("TF","Consensus","pvalue","anno")

foo[pvalue=="1"]
#foo[,log10pval:= -log10(pvalue)]
foo$log10pval<-do.call(rbind,strsplit(foo$pvalue,"[-]"))[,2]
foo[pvalue=="1",log10pval:=0]
foo$log10pval<-as.numeric(foo$log10pval)
foo$type.ndr<-do.call(rbind,strsplit(foo$anno,"[.]"))[,1]
foo$type.cor<-do.call(rbind,strsplit(foo$anno,"[.]"))[,2]
foo<-foo[type.cor%in% c("pos","neg")]
foo$stage<-do.call(rbind,strsplit(foo$anno,"[.]"))[,4]


kept.TF<-foo[,id:=paste0(TF,"_",Consensus)] %>%  .[log10pval>=10,unique(id)]
foo<-foo[id %in% kept.TF] %>% .[,id:=NULL]

### all may be used homer linked name consensus with tf
tf.anno.all<-fread("~/Desktop/NGS_postdoc/data/tfactivity/homer/20190806.homer.linked.NAME.consensus.motifid.tf.tsv") %>% setnames(c("TF","Consensus","id","motif.id","tf","extra.info"))
tf.anno.all[TF=="RAR:RXR(NR),DR5/ES-RAR-ChIP-Seq(GSE56893)/Homer"]
## the result from found consensus of tf may be written in different way, here we just link TF experiment name with tf symbol, so we don't need to consider the type of consensus
foo<-merge(foo,tf.anno.all[,.(TF,tf,extra.info)],by=c("TF")) %>% unique()
foo[extra.info=="unused"]


rna.tpm<-fread("~/Desktop/NGS_postdoc/data/RNA.expr.mat/20190519.tech.rep.collapsed.RNA.matrix.tpm.tsv")
gene.anno<-fread("~/Desktop/NGS_postdoc/data/DNA.data/ref/wy.gencode.gene.metadata.tsv")%>% .[,.(gene=ens_ID,symbol)]
sample.anno<-fread("~/Desktop/NGS_postdoc/data/DNA.data/20190524.updated.cell.metadata.tsv") %>% .[qc.rna==T]
unique(sample.anno$stage)
rna.tpm<-merge(rna.tpm,sample.anno[,.(cellid,stage)],by="cellid")
rna.tpm[,expr:=log2(tpm+1)]
rna.tpm.stage<-rna.tpm[,mean(tpm),.(stage,gene)]
rna.tpm.stage<-merge(rna.tpm.stage,gene.anno[,.(gene,symbol)],by="gene")
rna.tpm.stage[,symbol:=toupper(symbol)]
rna.tpm.stage[,unique(stage)]
tpm.ov.5<-rna.tpm.stage[V1>=5,unique(symbol)]

##  at least one stage tf expression >=5
foo<-foo[tf %in% tpm.ov.5]  
## for the duplicated tf , just keep the most significant one
test<-foo[,.N,tf]
test[N>8]

1: NR5A2 16

foo[tf=="NR5A2"]
foo<-foo[!TF=="Nr5a2(NR)/Pancreas-LRH1-ChIP-Seq(GSE34295)/Homer"]


### NOTE need to remove duplicated genes in DEG
deg<-fread("~/Desktop/NGS_postdoc/data/DEGs/20190610.degs/20190610.degs.TE.vs.ICM.deg.xls") %>% 
  .[,gene:=substring(V1,1,18)] %>% .[gene%in% gene.anno$gene] %>% 
  .[,.(symbol=toupper(substring(V1,20,100)),log2FoldChange,pvalue)]

foo[!tf %in% deg[,symbol],unique(tf)] #"NKX3-2" "SP5"    "THRB"  
#foo<-merge(foo,rna.tpm.stage[stage%in% "2cell",.(expr=log2(V1+1),tf=symbol)],by=c("tf"))
foo<-merge(foo,deg[,.(expr=log2FoldChange,tf=symbol,DEG.pvalue=pvalue)],by=c("tf"))
foo$type.ndr<-factor(foo$type.ndr,levels=c("promoter","distal"))#"zygote" ,"2cell", "4cell", "L4cell", "8cell","16cell","ICM" ,"TE"))
foo$stage<-factor(foo$stage,levels =rev( c("ICM","TE")))
foo$type.cor<-factor(foo$type.cor,levels=c("pos","neg"))#"zygote" ,"2cell", "4cell", "L4cell", "8cell","16cell","ICM" ,"TE"))
#foo$type.cor<-factor(foo$type.cor,levels=c("pos","neg","not"))#"zygote" ,"2cell", "4cell", "L4cell", "8cell","16cell","ICM" ,"TE"))

########

############################
foo[,max(log10pval)]
foo[,max(expr)]

to.plot<-copy(foo)#[ndrs=="allstage"]) 
to.plot[expr>2,expr:=2]
to.plot[expr< -2,expr:= -2]
to.plot[log10pval>100,log10pval:=100]
to.plot[,length(unique(tf))]

test<-dcast(to.plot,tf~type.ndr+stage+type.cor,value.var = "log10pval") %>% as.data.frame()%>% tibble::column_to_rownames("tf")
ord.tf<-pheatmap::pheatmap(test,clustering_distance_rows="correlation",cluster_cols = F)#,clustering_distance_rows="correlation")#,clustering_method = "ward.D2")

dev.off()
to.plot$tf<-factor(to.plot$tf,levels = rev(ord.tf$tree_row$labels[ord.tf$tree_row$order] ))
pdf("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/homer.TF/20190815.TF.enrichment.unmasked.sizegiven.pdf",width = 8,height=5)
ggplot(data = to.plot,aes(x=tf,y=stage))+
  geom_point(aes(size=log10pval,color=expr))+ 
  scale_color_gradient2(low="royalblue", midpoint = 0, mid = "grey",high = "orangered")+# rainbow(5)) + "grey90", "#",
  scale_size_continuous(breaks = c(0,10,25,50,100),range = c(0.2,4))+
  facet_grid(type.cor+type.ndr~.)+
  labs(y="",x="",title = "TF enrichment",color="Differential Expression\nlog2(fold change)",size="Motif Enrichment\n-log10(P value)")+
  theme_minimal()+
  theme(
    axis.text.y = element_text(size=10),
    #    axis.title.x = element_text(size=16),
    axis.title = element_text(size=14),
    panel.grid = element_blank(),
    #axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,size=14),
    axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,size=12),
    legend.title = element_text(size=12),
    legend.text = element_text(size=10),
    legend.position = "bottom",
    legend.box = "vertical",
    strip.text = element_text(size=10)
  )+ coord_fixed(ratio = 1)## coord_flip()
dev.off()

fwrite(foo[,.(TF= str_to_title(tf),Type=type.ndr, CREs=paste0(stage,".CRE.",type.cor,".sig"), `Motif Enrichment -log10(P value)`=log10pval,`Differential Expression log2(fold change)`=expr, `DEG P value`=DEG.pvalue)],
       "~/Desktop/Figure5g.tsv",sep="\t")

# save tf
fwrite(unique(foo[,.(tf,TF,Consensus)]),"~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/homer.TF/20190815.kept.homer.motif.tsv",sep="\t")

to.plot[log10pval<10,log10pval:=0]
to.plot[,length(unique(tf))]

test<-dcast(to.plot,tf~type.ndr+stage+type.cor,value.var = "log10pval") %>% as.data.frame()%>% tibble::column_to_rownames("tf")
ord.tf<-pheatmap::pheatmap(test,clustering_distance_rows="correlation",cluster_cols = F)#,clustering_distance_rows="correlation")#,clustering_method = "ward.D2")

dev.off()
to.plot$tf<-factor(to.plot$tf,levels = rev(ord.tf$tree_row$labels[ord.tf$tree_row$order] ))

pdf("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/homer.TF/20190815.TF.enrichment.unmasked.sizegiven.less10as0.pdf",width = 8,height=5)
ggplot(data = to.plot,aes(x=tf,y=stage))+
  geom_point(aes(size=log10pval,color=expr))+ 
  scale_color_gradient2(low="royalblue", midpoint = 0, mid = "grey",high = "orangered")+# rainbow(5)) + "grey90", "#",
  scale_size_continuous(breaks = c(0,10,25,50,100),range = c(0.2,4))+
  facet_grid(type.cor+type.ndr~.)+
  labs(y="",x="",title = "TF enrichment",color="Differential Expression\nlog2(fold change)",size="Motif Enrichment\n-log10(P value)")+
  theme_minimal()+
  theme(
    axis.text.y = element_text(size=10),
    #    axis.title.x = element_text(size=16),
    axis.title = element_text(size=14),
    panel.grid = element_blank(),
    #axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,size=14),
    axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,size=12),
    legend.title = element_text(size=12),
    legend.text = element_text(size=10),
    legend.position = "bottom",
    legend.box = "vertical",
    strip.text = element_text(size=10)
  )+ coord_fixed(ratio = 1)## coord_flip()
dev.off()

