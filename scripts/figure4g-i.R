library(data.table)
library(purrr)
library(furrr)
library(dplyr)
library(plyr)
### after run 20190714.sc.tf.activity.R
### link motif name to TF names
######################################################################################################
## tf activity in single cells
#figure 4f

zga.tf<-fread("~/Desktop/NGS_postdoc/data/zga/degs.cor/homer.tf/20190803.kept.homer.motif.tsv")  %>% setnames(c("TF","tf","Consensus"))

#zga.tf[,tmp.tf:=paste0(tf,"_",Consensus)]
tf.anno.all<-fread("~/Desktop/NGS_postdoc/data/tfactivity/homer/20190806.homer.linked.NAME.consensus.motifid.tf.tsv")
tf.anno.all<-tf.anno.all[tf%in% zga.tf$tf]
zga.tf
## there is duplicated tf with different consensus, after the motif enrichment homer result consensus may be written differ the original one
## so if use consensus for overlap should be very careful
tf.anno.all<-tf.anno.all[!motif.id=="motif_kt_56"]


files<-list.files("~/Desktop/NGS_postdoc/data/tfactivity/homer/sc.acc.with.bg",pattern = "20190805.TFBS.distal_NDR.ex100bp.bed.",full.names = T)

foo<-lapply(files,function(f)fread(f))%>% rbindlist() %>% .[,tf:=NULL]
foo<-merge(foo,tf.anno.all[,.(motif.id,TF,Consensus,extra.info)],by="motif.id")
foo[TF %in% zga.tf$TF, unique(TF)]

# overlap with sample and annotate stage and other informations, only use qc.acc.extra==T 
foo<-merge(foo,sample.anno[qc.rna==T &qc.dna==T &qc.acc.extra==T,.(cellid,stage,sample=DNA.info,gch_mean_rate,wcg_mean_rate)],by="sample")

to.plot<-copy(foo)#[ndrs=="allstage"]) 
to.plot[expr>2,expr:=2]
to.plot[expr< -2,expr:= -2]
to.plot[log10pval>100,log10pval:=100]
#to.plot[,log10pval:=round(log10pval)]
#to.plot[,log10pval:=log2(log10pval+1)]
to.plot[,length(unique(tf))]
#to.plot<-to.plot[type.cor=="pos"]

test<-dcast(to.plot,tf~type.ndr+type.cor,value.var = "log10pval") %>% as.data.frame()%>% tibble::column_to_rownames("tf")
ord.tf<-pheatmap::pheatmap(test,cluster_cols = F)#,clustering_distance_rows="correlation")#,clustering_method = "ward.D2")

dev.off()
#to.plot$tf<-factor(to.plot$tf,levels = rev(ord.tf$tree_row$labels[ord.tf$tree_row$order] ))


foo$TF<-factor(foo$TF,levels =rev (ord.tf$tree_row$labels[ord.tf$tree_row$order] ))
foo<-foo[stage %in% c("zygote","2cell")]
foo$stage<-factor(foo$stage,levels = c("zygote","2cell"))

ggplot(foo,aes(stage,zscore))+
  facet_grid(.~TF)+
  geom_violin(aes(color=stage))+
  geom_boxplot(width=0.1,outlier.color = NA)+
  geom_hline(yintercept = 0,linetype="dashed")+
  theme_classic()+
  theme(
    axis.text.x = element_blank()
  )

foo.dcast<-dcast(foo[stage%in%c("zygote","2cell")],TF~cellid,value.var = "zscore")%>% as.data.frame() %>% tibble::column_to_rownames("TF")
sample.anno<-fread("~/Desktop/NGS_postdoc/data/DNA.data/20190524.updated.cell.metadata.tsv")
sample.anno$stage<-factor(sample.anno$stage,levels = c("zygote","2cell","4cell","L4cell","8cell","16cell","ICM","TE"))
colanno<-sample.anno[cellid %in% colnames(foo.dcast),.(cellid,stage,DC.order)] %>% 
  setkey(stage,DC.order) %>% 
  .[,DC.order:=NULL] %>% as.data.frame() %>% tibble::column_to_rownames("cellid")
p1<-pheatmap::pheatmap(foo.dcast[,rownames(colanno)],annotation_col = colanno,#scale="row",
                       cluster_cols = F,cluster_rows = T,cellheight =10,cellwidth = 10,fontsize_row = 10,
                       show_colnames = F, show_rownames = T, border_color = "grey30",
                       breaks = seq(-20,20,by=4),
                       color = colorRampPalette(c("royalblue","black","orangered"))(length( seq(-20,20,by=4))),
                       annotation_colors=list(stage = c( zygote="#7FC97F",`2cell`= "#BEAED4")) #, `4cell`= "#FDC086",L4cell= "#808000",
                       #`8cell`= "#386CB0",`16cell`= "#650136", ICM= "#E64B35CC",TE= "#4DBBD5CC"
)
dev.off()

pdf("~/Desktop/NGS_postdoc/data/zga/degs.cor/homer.tf/20190807.zga.sc.tf.activity.distal.pdf",width=6,height = 5)
pheatmap::pheatmap(foo.dcast[p1$tree_row$labels[p1$tree_row$order] ,rownames(colanno)],annotation_col = colanno,#scale="row",
                   cluster_cols = F,cluster_rows = F,cellheight =10,cellwidth = 10,fontsize_row = 10,
                   show_colnames = F, show_rownames = T, border_color = "grey30",
                   breaks = seq(-20,20,by=4),
                   color = colorRampPalette(c("royalblue","black","orangered"))(length( seq(-20,20,by=4))), main = "TF activity (Distal)",
                   annotation_colors=list(stage = c( zygote="#7FC97F",`2cell`= "#BEAED4")) #, `4cell`= "#FDC086",L4cell= "#808000",
                   #`8cell`= "#386CB0",`16cell`= "#650136", ICM= "#E64B35CC",TE= "#4DBBD5CC"
)
dev.off()



################### plot tf expression
rna.tpm<-merge(rna.tpm,gene.anno,by="gene") %>% .[,symbol:=toupper(symbol)]

to.plot<-rna.tpm[symbol %in% rownames(foo.dcast) & cellid %in% colnames(foo.dcast)] %>% 
  dcast(symbol~cellid,value.var="expr") %>% 
  as.data.frame() %>% 
  tibble::column_to_rownames("symbol")
pdf("~/Desktop/NGS_postdoc/data/zga/degs.cor/homer.tf/20190807.zga.sc.tf.expression.pdf",width=6,height = 5)
pheatmap::pheatmap(to.plot[p1$tree_row$labels[p1$tree_row$order] ,rownames(colanno)],annotation_col = colanno,#scale="row",
                   cluster_cols = F,cluster_rows = F,cellheight =10,cellwidth = 10,fontsize_row = 10,
                   show_colnames = F, show_rownames = T, border_color = "grey30",
                   breaks = seq(0,10,by=1),
                   color = colorRampPalette(c("white","orangered"))(length( seq(0,10,by=1))), main = "TF expression",
                   annotation_colors=list(stage = c( zygote="#7FC97F",`2cell`= "#BEAED4")) #, `4cell`= "#FDC086",L4cell= "#808000",
                   #`8cell`= "#386CB0",`16cell`= "#650136", ICM= "#E64B35CC",TE= "#4DBBD5CC"
)
dev.off()



######################################################################################################################
## tf activity in single cells
files<-list.files("~/Desktop/NGS_postdoc/data/zga/degs.cor/homer.tf",recursive = T,pattern = "Results.txt",full.names = T)
basename(dirname(files))
foo<-lapply(files, function(f)fread(f,select = c(1,2,3)) %>% .[,anno:=basename(dirname(f)) %>% 
                                                                 sub("20190802.accrna.correlations.","",.)%>%
                                                                 sub(".tsv_motif_unmask_sizeGiven","",.)])%>% rbindlist()

colnames(foo)<-c("TF","Consensus","pvalue","anno")

#foo$tf<-do.call(rbind,strsplit(foo$TF,"[/]"))[,1]
#foo$tf<-do.call(rbind,strsplit(foo$tf,"[()]"))[,1]
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

tf.anno.all<-fread("~/Desktop/NGS_postdoc/data/tfactivity/homer/20190806.homer.linked.NAME.consensus.motifid.tf.tsv")
zga.tf<-fread("~/Desktop/NGS_postdoc/data/zga/degs.cor/homer.tf/20190803.kept.homer.motif.tsv")  %>% setnames(c("TF","tf","Consensus"))

## there is duplicated tf with different consensus, after the motif enrichment homer result consensus may be written differ the original one
## so if use consensus for overlap should be very careful
## for duplicated tf(have different consensus)，keep most sig hit
## homer consensus may differ from ref

figure 4g

foo[,len:=nchar(Consensus)]
tf.anno.all[,len:=nchar(Consensus)]
foo<-merge(foo[,.(tf=TF,Consensus,len,log10pval,anno)],tf.anno.all,by=c("tf","len"))
tf.anno.all[,length(unique(TF))]  #354
kept.tf<-foo[,.SD[which.max(log10pval)],TF]
kept.tf[,length(unique(TF))] #352
kept.tf[tf%in%zga.tf$tf,length(unique(TF))] #22
zga.tf[!tf%in% kept.tf$tf]

## load distal sc tf acc 
files<-list.files("~/Desktop/NGS_postdoc/data/tfactivity/homer/sc.acc.with.bg",pattern = "20190805.TFBS.distal_NDR.ex100bp.bed.",full.names = T)
foo<-lapply(files,function(f)fread(f))%>% rbindlist() %>% .[,tf:=NULL]
## transfer motif.id to TF, rm duplicated tf symbol, only keep most sig hit
foo<-merge(foo,kept.tf[,.(motif.id,TF)],by="motif.id")
foo[TF %in% zga.tf$TF, unique(TF)]

# overlap with sample and annotate stage and other informations, only use qc.acc.extra==T 
foo<-merge(foo,sample.anno[qc.rna==T &qc.dna==T &qc.acc.extra==T,.(cellid,stage,sample=DNA.info,gch_mean_rate,wcg_mean_rate)],by="sample")

#foo$TF<-factor(foo$TF,levels =rev (ord.tf$tree_row$labels[ord.tf$tree_row$order] ))
foo$stage<-factor(foo$stage,levels = c("zygote","2cell","4cell","L4cell","8cell","16cell","ICM","TE"))
kept.tf[motif.id=="motif_kt_308"]
kept.tf[TF=="RARA"]
kept.tf[TF=="RARG"]
kept.tf[TF=="RARG"]

kept.tf[motif.id=="motif_kt_308"]

tmp<-foo[ weight>=50 &stage%in% c("zygote","2cell")]# &cellid %in% sample.anno[acc.type=="hiAccRNA_loAcc",cellid]]
tmp<-tmp[,t.test(zscore~stage),.(TF,motif.id)]
tmp[,"stage":=c("zygote","2cell"),TF]
tmp<-dcast(tmp,TF+motif.id+p.value~stage,value.var = "estimate") %>% setkey(p.value)
tmp[p.value<0.05]
tmp[,padj:=p.adjust(p.value, method="fdr")]
tmp[padj<0.1 ]
tmp[TF=="RARA"]

tmp[TF %in%zga.tf$TF,plot.color:="enriched"]
tmp[TF %in%zga.tf$TF &padj<0.1,plot.color:="enriched_sig"]
#tmp[!TF %in%zga.tf$TF,plot.color:="not.enriched"]
#tmp[!TF %in%zga.tf$TF&padj<0.1,plot.color:="not.enriched_sig"]
setorder(tmp,-padj)

pdf("~/Desktop/NGS_postdoc/data/zga/degs.cor/20190820.sc.TF.activity.diff.pdf",width=5,height=5)
xlen<-tmp[TF %in% zga.tf$TF,max(zygote)-min(zygote)]
ylen<-tmp[TF %in% zga.tf$TF,max(`2cell`)-min(`2cell`)]
tmp[TF %in% zga.tf$TF,.SD[which.max(`2cell`)]]
ggplot(tmp[TF %in% zga.tf$TF],aes(zygote,`2cell`,color=plot.color))+
  geom_point(alpha=0.8,size=1)+
  scale_color_manual(values = c(enriched="black",enriched_sig="orangered"))+
  geom_abline(slope = 1,alpha=0.5)+
  labs(x="TF activity in zygote",y="TF activity in 2cell")+
  geom_hline(yintercept = 0,linetype="dashed",alpha=0.5)+
  geom_vline(xintercept = 0,linetype="dashed",alpha=0.5)+
  theme_bw()+
  ggrepel::geom_text_repel(data=tmp[TF %in%zga.tf$TF &padj<0.1 &p.value<0.05],aes(x=zygote, y=`2cell`, label=TF), vjust=-0.0, hjust=-0.3,size=3,color="red")+ 
#  ggrepel::geom_text_repel(data=tmp[TF %in% c("GSC","OTX2")],aes(x=zygote, y=`2cell`, label=TF), vjust=-0.5, hjust=-0.5,size=3,color="black") +
  theme(
    panel.grid = element_blank()
  )+coord_fixed(ratio=xlen/ylen)
dev.off()

diff.tf.activity<-tmp[plot.color=="enriched_sig"]




### tf expression
degs<-fread("~/Desktop/NGS_postdoc/data/DEGs/20190610.zga.gene.fc.4.padj0.01.basemean10.tsv")
degs<-degs[,symbol:=toupper(symbol)]
degs[symbol=="OTX2"]

rna.tpm<-fread("~/Desktop/NGS_postdoc/data/RNA.expr.mat/20190519.tech.rep.collapsed.RNA.matrix.tpm.tsv")
rna.tpm<-merge(rna.tpm,sample.anno[qc.rna==T& !sex=="GG" & stage %in% c("zygote","2cell"),.(cellid,stage)],by="cellid")
rna.tpm[,expr:=log2(tpm+1)]
rna.tpm<-merge(rna.tpm,gene.anno[,.(gene,symbol)],by="gene")
rna.tpm[,symbol:=toupper(symbol)]
## cal diff expr between zygote and 2cell 
tmp<-rna.tpm[symbol%in%zga.tf$TF]
filter.low.expr.genes<-tmp[,mean(expr),.(symbol,stage)]%>% .[V1>1,unique(symbol)]
tmp<-tmp[symbol %in% filter.low.expr.genes,t.test(expr~stage),.(symbol)]
tmp[symbol=="GSC"]
tmp[,"stage":=c("2cell","zygote"),symbol]
tmp<-dcast(tmp[!is.nan(p.value)],symbol+p.value~stage,value.var = "estimate") %>% setkey(p.value)
tmp[p.value<0.05]
tmp[,padj:=p.adjust(p.value, method="fdr")]
tmp[padj<0.01 ]


tmp[symbol %in%zga.tf$TF,plot.color:="enriched"]
tmp[symbol %in%zga.tf$TF &padj<0.01,plot.color:="enriched_sig"]
#tmp[!symbol %in%zga.tf$TF,plot.color:="not.enriched"]
#tmp[!symbol %in%zga.tf$TF&padj<0.01,plot.color:="not.enriched_sig"]
setorder(tmp,-padj)

tmp[symbol=="OTX2"]

pdf("~/Desktop/NGS_postdoc/data/zga/degs.cor/20190820.sc.TF.expr.diff.pdf",width=5,height=5)
xlen<-tmp[,max(zygote)-min(zygote)]
ylen<-tmp[,max(`2cell`)-min(`2cell`)]
tmp[,.SD[which.max(`2cell`)]]
ggplot(tmp,aes(zygote,`2cell`,color=plot.color))+
  geom_point(alpha=0.8,size=1)+
  scale_color_manual(values = c(enriched="black",enriched_sig="orangered"))+
  geom_abline(slope = 1,alpha=0.5)+
  labs(x="zygote log2(TPM+1)",y="2cell log2(TPM+1)")+
  #geom_hline(yintercept = 0,linetype="dashed",alpha=0.5)+
  #geom_vline(xintercept = 0,linetype="dashed",alpha=0.5)+
  theme_bw()+
  ggrepel::geom_text_repel(data=tmp[plot.color=="enriched_sig"],aes(x=zygote, y=`2cell`, label=symbol), vjust=-0.0, hjust=-0.3,size=3,color="red")+ 
#  ggrepel::geom_text_repel(data=tmp[symbol %in% c("GSC","OTX2")],aes(x=zygote, y=`2cell`, label=symbol), vjust=-0.5, hjust=-0.5,size=3,color="black") +
  theme(
    panel.grid = element_blank()
  )+coord_fixed(ratio=xlen/ylen)
dev.off()

diff.tf.expr<-tmp[plot.color=="enriched_sig"]


### merge expr and tf activity show sc 
to.plot<-merge(rna.tpm,foo[,.(cellid,symbol=TF,zscore)],by=c("cellid","symbol"))

kept.genes<-to.plot[,mean(tpm),.(symbol,stage)] %>% .[V1>1,unique(symbol)]
to.plot[,cor.test(expr,zscore),symbol]

cor_samples <- to.plot[, .(V1 = unlist(cor.test(expr, zscore, alternative = "two.sided", method = "pearson")[c("estimate", "statistic", "p.value")])), by = c("symbol")] 

cor_samples <- cor_samples %>% .[, para := rep(c("r","t","p"), .N/3)] %>% 
  data.table::dcast(symbol ~ para, value.var = "V1") %>%
  .[, c("padj_fdr", "padj_bonf") := list(p.adjust(p, method="fdr"), p.adjust(p, method="bonferroni"))] %>%
  .[, c("log_padj_fdr","log_padj_bonf") := list(-log10(padj_fdr), -log10(padj_bonf))] %>%
  .[, sig := padj_fdr <0.1] %>% 
  setorder(padj_fdr)
cor_samples[sig==T &symbol%in% zga.tf$TF]
cor_samples[symbol=="RARG"]

setorder(cor_samples,-padj_fdr)

#figure 4h
#### correlation plot
pdf("~/Desktop/NGS_postdoc/data/zga/degs.cor/20190820.sc.TF.expr.activity.pearson.cor.pdf",width=8,height=5)
xlen<-cor_samples[!is.na(r),max(r)-min(r)]
ylen<-cor_samples[!is.na(r),max(-log10(p))-min(-log10(p))]
ggplot(cor_samples[!is.na(r)],aes(r,-log10(p)))+
  geom_point(aes(color= (sig==T & symbol %in% zga.tf$TF)),alpha=0.8,size=0.5)+
  scale_color_manual(values = c("FALSE"="grey","TRUE"="red"))+
  ggrepel::geom_text_repel(data=cor_samples[sig==T &symbol %in% zga.tf$TF],aes(x=r, y=-log10(p), label=symbol), vjust=-0.0, hjust=-0.3,size=5,color="red")+
#  geom_vline(xintercept = 0)+
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )+coord_fixed(ratio=xlen/ylen)
dev.off()


figure 4i
diff.tf.activity
diff.tf.expr[symbol %in% diff.tf.activity$TF]
symbol%in%diff.tf.activity$TF |symbol%in% diff.tf.expr$symbol

for (i in diff.tf.expr[symbol %in% diff.tf.activity$TF,symbol] ) {
  pdf(sprintf("~/Desktop/NGS_postdoc/data/zga/degs.cor/20190821.sc.point.diff.acc_diff.expr_%s.pdf",i),width=4,height=3)
  x.len= to.plot[symbol==i,max(expr)]- to.plot[symbol==i,min(expr)]
  y.len= to.plot[symbol==i,max(zscore)]- to.plot[symbol==i,min(zscore)]
  r= cor_samples[symbol==i,round(r*100)/100]
  p1<-ggplot(to.plot[symbol==i],aes(expr,zscore))+
    geom_point(size=1,alpha=0.8,aes(color=stage))+
    scale_color_manual(values = c(  zygote="#7FC97F",`2cell`= "#BEAED4"))+
    labs(x="TF expression log2(TPM+1)",y="TF acctivity",title =paste0(i ," r=",r) )+
    geom_smooth(method = lm,se=F)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(size=16),
      axis.text = element_text(size=12)
    )+coord_fixed(ratio=x.len/y.len)
  print(p1)
  dev.off()
}


for (i in diff.tf.expr[!symbol %in% diff.tf.activity$TF,symbol] ) {
  pdf(sprintf("~/Desktop/NGS_postdoc/data/zga/degs.cor/20190821.sc.point.not.diff.acc_diff.expr_%s.pdf",i),width=4,height=3)
  x.len= to.plot[symbol==i,max(expr)]- to.plot[symbol==i,min(expr)]
  y.len= to.plot[symbol==i,max(zscore)]- to.plot[symbol==i,min(zscore)]
  r= cor_samples[symbol==i,round(r*100)/100]
  p1<-ggplot(to.plot[symbol==i],aes(expr,zscore))+
    geom_point(size=1,alpha=0.8,aes(color=stage))+
    scale_color_manual(values = c(  zygote="#7FC97F",`2cell`= "#BEAED4"))+
    labs(x="TF expression log2(TPM+1)",y="TF acctivity",title =paste0(i ," r=",r) )+
    geom_smooth(method = lm,se=F)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(size=16),
      axis.text = element_text(size=12)
    )+coord_fixed(ratio=x.len/y.len)
  print(p1)
  dev.off()
}

for (i in diff.tf.activity[!TF %in% diff.tf.expr$symbol,TF] ) {
  pdf(sprintf("~/Desktop/NGS_postdoc/data/zga/degs.cor/20190821.sc.point.diff.acc_not.diff.expr_%s.pdf",i),width=4,height=3)
  x.len= to.plot[symbol==i,max(expr)]- to.plot[symbol==i,min(expr)]
  y.len= to.plot[symbol==i,max(zscore)]- to.plot[symbol==i,min(zscore)]
  r= cor_samples[symbol==i,round(r*100)/100]
  p1<-ggplot(to.plot[symbol==i],aes(expr,zscore))+
    geom_point(size=1,alpha=0.8,aes(color=stage))+
    scale_color_manual(values = c(  zygote="#7FC97F",`2cell`= "#BEAED4"))+
    labs(x="TF expression log2(TPM+1)",y="TF acctivity",title =paste0(i ," r=",r) )+
    geom_smooth(method = lm,se=F)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(size=16),
      axis.text = element_text(size=12)
    )+coord_fixed(ratio=x.len/y.len)
  print(p1)
  dev.off()
}


