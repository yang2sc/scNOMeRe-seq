# figure 4e
# tf enrichment
library(data.table)
library(purrr)
library(furrr)
library(dplyr)
library(plyr)


#run homer in server, size given; plot in local
files<-list.files("~/Desktop/NGS_postdoc/data/zga/degs.cor/homer.tf",recursive = T,pattern = "Results.txt",full.names = T)
basename(dirname(files))
foo<-lapply(files, function(f)fread(f,select = c(1,2,3)) %>% .[,anno:=basename(dirname(f)) %>% 
                                                                 sub("20190802.accrna.correlations.","",.)%>%
                                                                 sub(".tsv_motif_unmask_sizeGiven","",.)])%>% rbindlist()

colnames(foo)<-c("TF","Consensus","pvalue","anno")

str(foo)
foo[pvalue=="1"]
#foo[,log10pval:= -log10(pvalue)]
foo$log10pval<-do.call(rbind,strsplit(foo$pvalue,"[-]"))[,2]
foo[pvalue=="1",log10pval:=0]
foo$log10pval<-as.numeric(foo$log10pval)
#foo$stage<-do.call(rbind,strsplit(foo$anno,"[.]"))[,2]
foo$type.ndr<-do.call(rbind,strsplit(foo$anno,"[.]"))[,1]
foo[,type.cor:=anno %>% sub("distal.|promoter.","",.)]
foo<-foo[type.cor%in% c("pos.sig","neg.sig")]


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
tpm.ov.5<-rna.tpm.stage[stage%in% c("zygote","2cell") &V1>=5,unique(symbol)]

## in zygote 2cell stage at least one stage tf expression >=5
foo<-foo[tf %in% tpm.ov.5]  
## for the duplicated tf , just keep the most significant one
test<-foo[,.N,tf]
test[N>4]

foo[tf=="RARA"]
foo<-foo[!TF=="RARa(NR)/K562-RARa-ChIP-Seq(Encode)/Homer"]
foo[tf=="NR5A2"]
foo<-foo[!TF=="Nr5a2(NR)/Pancreas-LRH1-ChIP-Seq(GSE34295)/Homer"]


### NOTE need to remove duplicated genes in DEG
deg<-fread("~/Desktop/NGS_postdoc/data/DEGs/20190610.degs/20190610.degs.2cell.vs.zygote.deg.xls") %>% 
  .[,gene:=substring(V1,1,18)] %>% .[gene%in% gene.anno$gene] %>% 
  .[,.(symbol=toupper(substring(V1,20,100)),log2FoldChange,pvalue)]

foo[!tf %in% deg[,symbol],unique(tf)] #
foo<-merge(foo,deg[,.(expr=log2FoldChange,tf=symbol,DEG.pvalue=pvalue)],by=c("tf"))
foo$type.ndr<-factor(foo$type.ndr,levels=c("promoter","distal"))
foo$type.cor<-factor(foo$type.cor,levels=c("pos.sig","neg.sig"))
fwrite(foo[,.(TF= str_to_title(tf),Type=type.ndr, CREs=type.cor, `Motif Enrichment -log10(P value)`=log10pval,`Differential Expression log2(fold change)`=expr, `DEG P value`=DEG.pvalue)],
       "~/Desktop/Figure4e.tsv",sep="\t")

####################################
foo[,max(log10pval)]
foo[,max(expr)]

to.plot<-copy(foo)#[ndrs=="allstage"]) 
to.plot[expr>2,expr:=2]
to.plot[expr< -2,expr:= -2]
to.plot[log10pval>100,log10pval:=100]
to.plot[,length(unique(tf))]
to.plot$type.cor<-factor(to.plot$type.cor,levels = c("neg.sig","pos.sig"))

test<-dcast(to.plot,tf~type.ndr+type.cor,value.var = "log10pval") %>% as.data.frame()%>% tibble::column_to_rownames("tf")
ord.tf<-pheatmap::pheatmap(test,cluster_cols = F)

dev.off()
to.plot$tf<-factor(to.plot$tf,levels = rev(ord.tf$tree_row$labels[ord.tf$tree_row$order] ))

pdf("~/Desktop/NGS_postdoc/data/zga/degs.cor/homer.tf/20190803.TF.enrichment.unmasked.sizegiven.pdf",width = 5,height=7)
ggplot(data = to.plot, mapping = aes(x=tf,y=type.cor))+
  geom_point(aes(size=log10pval,color=expr))+ 
  scale_color_gradient2(low="royalblue", midpoint = 0, mid = "grey",high = "orangered")+
  scale_size_continuous(breaks = c(0,10,25,50,100),range = c(0.2,4))+
  facet_grid(type.ndr~.)+
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


