# sup figure10b
## tf activity in single cells
files<-list.files("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/homer.TF",recursive = T,pattern = "Results.txt",full.names = T)
basename(dirname(files))
foo<-lapply(files, function(f)fread(f,select = c(1,2,3)) %>% .[,anno:=basename(dirname(f)) %>% 
                                                                 sub("20190810.accrna.correlations.","",.)%>%
                                                                 sub(".high.tsv_motif_unmask_sizeGiven","",.)])%>% rbindlist()
colnames(foo)<-c("TF","Consensus","pvalue","anno")

#foo$tf<-do.call(rbind,strsplit(foo$TF,"[/]"))[,1]
#foo$tf<-do.call(rbind,strsplit(foo$tf,"[()]"))[,1]
str(foo)
foo[pvalue=="1"]
#foo[,log10pval:= -log10(pvalue)]
foo$log10pval<-do.call(rbind,strsplit(foo$pvalue,"[-]"))[,2]
foo[pvalue=="1",log10pval:=0]
foo$log10pval<-as.numeric(foo$log10pval)
foo$type.cor<-do.call(rbind,strsplit(foo$anno,"[.]"))[,2]
foo<-foo[type.cor%in% c("pos","neg")]

tf.anno.all<-fread("~/Desktop/NGS_postdoc/data/tfactivity/homer/20190806.homer.linked.NAME.consensus.motifid.tf.tsv")
icm.te.tf<-fread("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/homer.TF/20190815.kept.homer.motif.tsv")  %>% setnames(c("TF","tf","Consensus"))

## there is duplicated tf with different consensus, after the motif enrichment homer result consensus may be written differ the original one
## so if use consensus for overlap should be very careful

foo[,len:=nchar(Consensus)]
tf.anno.all[,len:=nchar(Consensus)]
foo<-merge(foo[,.(tf=TF,Consensus,len,log10pval,anno)],tf.anno.all,by=c("tf","len"))
tf.anno.all[,length(unique(TF))]  #354
kept.tf<-foo[,.SD[which.max(log10pval)],TF]
kept.tf[,length(unique(TF))] #352
kept.tf[tf%in%icm.te.tf$tf,length(unique(TF))] #33
icm.te.tf[!tf%in% kept.tf$tf]

## load distal sc tf acc 
files<-list.files("~/Desktop/NGS_postdoc/data/tfactivity/homer/sc.acc.with.bg",pattern = "20190805.TFBS.distal_NDR.ex100bp.bed.",full.names = T)
foo<-lapply(files,function(f)fread(f))%>% rbindlist() %>% .[,tf:=NULL]
foo<-merge(foo,kept.tf[,.(motif.id,TF)],by="motif.id")
foo[TF %in% icm.te.tf$TF, unique(TF)]

# overlap with sample and annotate stage and other informations, only use qc.acc.extra==T 
foo<-merge(foo,sample.anno[qc.rna==T &qc.dna==T &qc.acc.extra==T,.(cellid,stage,sample=DNA.info,gch_mean_rate,wcg_mean_rate)],by="sample")

#foo$TF<-factor(foo$TF,levels =rev (ord.tf$tree_row$labels[ord.tf$tree_row$order] ))
foo$stage<-factor(foo$stage,levels = c("zygote","2cell","4cell","L4cell","8cell","16cell","ICM","TE"))
kept.tf[motif.id=="motif_kt_308"]
kept.tf[TF=="RARA"]
kept.tf[TF=="RARG"]
kept.tf[TF=="RARG"]

kept.tf[motif.id=="motif_kt_308"]



foo.dcast<-dcast(foo[TF %in% icm.te.tf$TF],TF~cellid,value.var = "zscore")%>% as.data.frame() %>% tibble::column_to_rownames("TF")
sample.anno<-fread("~/Desktop/NGS_postdoc/data/DNA.data/20190524.updated.cell.metadata.tsv")
sample.anno$stage<-factor(sample.anno$stage,levels = c("zygote","2cell","4cell","L4cell","8cell","16cell","ICM","TE"))
colanno<-sample.anno[cellid %in% colnames(foo.dcast),.(cellid,stage,DC.order)] %>% 
  setkey(stage,DC.order) %>% 
  .[,DC.order:=NULL] %>% as.data.frame() %>% tibble::column_to_rownames("cellid")

p1<-pheatmap::pheatmap(foo.dcast[,rownames(colanno)],annotation_col = colanno,#scale="row",
                       cluster_cols = F,cluster_rows = T,cellheight =10,cellwidth = 3,fontsize_row = 10,
                       show_colnames = F, show_rownames = T, border_color = NA,
                       breaks = seq(-20,20,by=2),
                       color = colorRampPalette(c("royalblue","white","orangered"))(length( seq(-20,20,by=2))),
                       annotation_colors=list(stage = c( zygote="#7FC97F",`2cell`= "#BEAED4" , `4cell`= "#FDC086",L4cell= "#808000",
                                                         `8cell`= "#386CB0",`16cell`= "#650136", ICM= "#E64B35CC",TE= "#4DBBD5CC"))
)
dev.off()
#  #(ord.tf$tree_row$labels[ord.tf$tree_row$order] )
pdf("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/homer.TF/20190816.icm.te.sc.tf.activity.distal.pdf",width=6,height = 5)
pheatmap::pheatmap(foo.dcast[p1$tree_row$labels[p1$tree_row$order],rownames(colanno)],annotation_col = colanno,#scale="row",
                   cluster_cols = F,cluster_rows = F,cellheight =10,cellwidth = 2,fontsize_row = 10,
                   show_colnames = F, show_rownames = T, border_color = "grey30",
                   breaks = seq(-20,20,by=2),
                   color = colorRampPalette(c("lightblue","black","orangered"))(length( seq(-20,20,by=2))), main = "TF activity (Distal)",
                   annotation_colors=list(stage = c( zygote="#7FC97F",`2cell`= "#BEAED4", `4cell`= "#FDC086",L4cell= "#808000",
                                                     `8cell`= "#386CB0",`16cell`= "#650136", ICM= "#E64B35CC",TE= "#4DBBD5CC")) 
)
dev.off()


################### plot tf expression
rna.tpm<-fread("~/Desktop/NGS_postdoc/data/RNA.expr.mat/20190519.tech.rep.collapsed.RNA.matrix.tpm.tsv")
gene.anno<-fread("~/Desktop/NGS_postdoc/data/DNA.data/ref/wy.gencode.gene.metadata.tsv")%>% .[,.(gene=ens_ID,symbol)]
rna.tpm<-merge(rna.tpm,sample.anno[qc.rna==T, .(cellid,stage,DC.order)],by="cellid") %>% merge(gene.anno,by="gene")
rna.tpm[,expr:=log2(tpm+1)]

rna.tpm<-rna.tpm[,symbol:=toupper(symbol)]

to.plot<-rna.tpm[symbol %in% rownames(foo.dcast) & cellid %in% colnames(foo.dcast)] %>% 
  dcast(symbol~cellid,value.var="expr") %>% 
  as.data.frame() %>% 
  tibble::column_to_rownames("symbol")
pdf("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/homer.TF/20190816.icm.te.sc.tf.expression.pdf",width=6,height = 5)
pheatmap::pheatmap(to.plot[p1$tree_row$labels[p1$tree_row$order] ,rownames(colanno)],annotation_col = colanno,#scale="row",
                   cluster_cols = F,cluster_rows = F,cellheight =10,cellwidth = 2,fontsize_row = 10,
                   show_colnames = F, show_rownames = T, border_color = NA,
                   breaks = seq(0,10,by=0.5),
                   color = colorRampPalette(c("black","yellow"))(length( seq(0,10,by=0.5))), main = "TF expression",
                   annotation_colors=list(stage = c( zygote="#7FC97F",`2cell`= "#BEAED4", `4cell`= "#FDC086",L4cell= "#808000",
                                                     `8cell`= "#386CB0",`16cell`= "#650136", ICM= "#E64B35CC",TE= "#4DBBD5CC"))
)
dev.off()




