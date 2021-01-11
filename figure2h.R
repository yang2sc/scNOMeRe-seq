# figure 2h

library(data.table)
## rna expression
sample.anno<-fread("~/Desktop/NGS_postdoc/data/DNA.data/20190826.updated.sample.metadata.tsv")
rna<-fread("~/Desktop/NGS_postdoc/data/RNA.expr.mat/20190519.tech.rep.collapsed.RNA.matrix.tpm.tsv") %>% merge(sample.anno[qc.rna==T &!sex=="GG",.(cellid,stage,embryo_id)],by="cellid")
gene.anno<-fread("~/Desktop/NGS_postdoc/data/DNA.data/ref/wy.gencode.gene.metadata.tsv") %>% .[,.(gene=ens_ID,symbol)]
rna<-merge(rna,gene.anno,by="gene")
rna[,expr:=log2(tpm+1)]

kept.gene1<-rna[,.(mean(expr)),.(stage,symbol)] %>% .[V1>1,unique(symbol)] #18525

#  cv >0.5
#kept.gene3<-rna[,sd(expr)/mean(expr),symbol] %>% .[V1>0.5,unique(symbol)] 
#rna[symbol %in% kept.gene1 & symbol %in% kept.gene3,length(unique(symbol))] #  10704
#correlation.test<-rna[symbol %in% kept.gene1 & symbol %in% kept.gene3] 

correlation.test<-rna [symbol %in% kept.gene1 ]

correlation.test<-dcast(correlation.test,symbol~cellid,value.var = "expr") %>% as.data.frame() %>% tibble::column_to_rownames("symbol")
correlation.test<-cor(correlation.test,method = "spearman")
#row.anno<-sample.anno[cellid %in% rownames(correlation.test) &qc.rna==T,.(cellid,DC.order,stage)]%>% setorder(DC.order) %>% as.data.frame() %>% tibble::column_to_rownames("cellid")
#pheatmap::pheatmap(correlation.test[rownames(row.anno),rownames(row.anno)],cluster_rows = F,cluster_cols = F,annotation_row = row.anno)

to.plot<-flattenCorrMatrix(correlation.test)%>% as.data.table %>% setnames(c("cell1","cell2","cor")) 
to.plot[,length(unique(cell2))]

to.plot<-merge(to.plot,sample.anno[,.(cell1=cellid,stage1=stage,embryo_id1=embryo_id,m_id1=mother_id,g_id1=g_id)],by="cell1", all.x=T)
to.plot<-merge(to.plot,sample.anno[,.(cell2=cellid,stage2=stage,embryo_id2=embryo_id,m_id2=mother_id,g_id2=g_id)],by="cell2", all.x=T)
unique(to.plot$embryo_id1)
to.plot[embryo_id1==embryo_id2 & embryo_id2=="X4cell2"]

#to.plot[stage1 == "ICM", stage1:= "blastocyst"]
#to.plot[stage1 == "TE", stage1:= "blastocyst"]
#to.plot[stage2 == "ICM", stage2:= "blastocyst"]
#to.plot[stage2 == "TE", stage2:= "blastocyst"]
unique(to.plot$stage1)
to.plot$stage1<-factor(to.plot$stage1,levels=c("zygote","2cell","4cell","L4cell","8cell","16cell", "ICM","TE")) #"ICM","TE"

to.plot[stage1==stage2, group:=ifelse(embryo_id1==embryo_id2,"intraembryo","interembryo")]
to.plot<-to.plot[stage1==stage2]
to.plot$group<-factor(to.plot$group,levels = c("intraembryo","interembryo"))


pdf("~/Desktop/NGS_postdoc/data/lineage_tracing/20190827.RNA.tpm.embryo.compare.box.pdf",width=5,height=3)
ggplot(to.plot[stage1==stage2],aes(stage1,cor,color=group))+
  geom_boxplot( outlier.size = 0.2,position = position_dodge(preserve = "single",width=0.9))+
  #geom_jitter()+
  ylim(0.7,1)+
  scale_color_manual(values=c("#DF8F44CC","#00A1D5CC"))+
  theme_bw()+
  labs(x="",y="Pairewised Spearman correaltions",title = "Transcriptome")+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=45,hjust=1),
    plot.title = element_text(hjust=0.5)
  )#+coord_fixed(ratio=4/0.7)
dev.off()

to.plot[,length(unique(cell1))]

to.plot[group=="intraembryo", mother_type:=ifelse( m_id1==m_id2,"same mother","different mother")]
to.plot$mother_type<-factor(to.plot$mother_type,levels = c("same mother","different mother"))
pdf("~/Desktop/NGS_postdoc/data/lineage_tracing/20190827.RNA.tpm.mother.compare.box.pdf",width=5,height=3)
ggplot(to.plot[group=="intraembryo" &stage1%in% c("2cell","4cell","L4cell","8cell")],aes(stage1,cor,color=mother_type))+
  geom_boxplot( outlier.size = 0.2,position = position_dodge(preserve = "single",width=0.9))+
  #geom_jitter()+
  ylim(0.7,1)+
  theme_bw()+
  scale_color_manual(values=c("#DC0000CC","#008280CC"))+
  labs(x="",y="Pairewised Spearman correaltions",title = "Transcriptome")+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=45,hjust=1),
    plot.title = element_text(hjust=0.5)
  )#+coord_fixed(ratio=4/0.7)
dev.off()

to.plot[group=="intraembryo", grand_type:=ifelse( g_id1==g_id2,"same grand","different grand")]
to.plot$grand_type<-factor(to.plot$grand_type,levels = c("same grand","different grand"))
pdf("~/Desktop/NGS_postdoc/data/lineage_tracing/20190827.RNA.tpm.grand.compare.box.pdf",width=5,height=3)
ggplot(to.plot[group=="intraembryo" &stage1%in% c("4cell","L4cell","8cell")],aes(stage1,cor,color=grand_type))+
  geom_boxplot( outlier.size = 0.2,position = position_dodge(preserve = "single",width=0.9))+
  #geom_jitter()+
  ylim(0.7,1)+
  scale_color_startrek()+
  theme_bw()+
  labs(x="",y="Pairewised Spearman correaltions",title = "Transcriptome")+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=45,hjust=1),
    plot.title = element_text(hjust=0.5)
  )#+coord_fixed(ratio=4/0.7)
dev.off()

to.plot<-rbind(to.plot[,.(Type=group,Stage=stage1,cor,Data="Expr")],
               to.plot[group=="intraembryo" &stage1%in% c("2cell","4cell","L4cell","8cell"),.(Type=mother_type,Stage=stage1,cor,Data="Expr")],
               to.plot[group=="intraembryo" &stage1%in% c("4cell","L4cell","8cell"),.(Type=grand_type,Stage=stage1,cor,Data="Expr")]
)

mypal = pal_jama("default", alpha = 0.8)(10)
mypal
library("scales")
show_col(mypal)
####    met 5kb bin 
## 20190717.clustering.R
#acc.cor.dcast<-read.table("~/Desktop/NGS_postdoc/data/clustering/20190716.acc.ndr.spearman.correlations.tsv")
#rna.cor.dcast<-read.table("~/Desktop/NGS_postdoc/data/clustering/20190716.expr.tpm.spearman.correlations.tsv")

correlation.test<-read.table("~/Desktop/NGS_postdoc/data/clustering/20190716.met.5kb.spearman.correlations.tsv")

to.plot<-flattenCorrMatrix(correlation.test)%>% as.data.table %>% setnames(c("sample1","sample2","cor")) 
to.plot[,length(unique(sample1))]
#### remove pg cells
to.plot<-merge(to.plot,sample.anno[!sex=="GG",.(sample1=DNA.info,cell1=cellid,stage1=stage,embryo_id1=embryo_id,m_id1=mother_id,g_id1=g_id)],by="sample1")
to.plot<-merge(to.plot,sample.anno[!sex=="GG",.(sample2=DNA.info,cell2=cellid,stage2=stage,embryo_id2=embryo_id,m_id2=mother_id,g_id2=g_id)],by="sample2")
sample.anno[cellid%in% to.plot$cell2,.N,sex]

unique(to.plot$embryo_id1)
to.plot[embryo_id1==embryo_id2 & embryo_id2=="X4cell2"]

#to.plot[stage1 == "ICM", stage1:= "blastocyst"]
#to.plot[stage1 == "TE", stage1:= "blastocyst"]
#to.plot[stage2 == "ICM", stage2:= "blastocyst"]
#to.plot[stage2 == "TE", stage2:= "blastocyst"]
unique(to.plot$stage1)
to.plot$stage1<-factor(to.plot$stage1,levels=c("zygote","2cell","4cell","L4cell","8cell","16cell", "ICM","TE" )) #"blastocyst")) #

to.plot[stage1==stage2, group:=ifelse(embryo_id1==embryo_id2,"intraembryo","interembryo")]
to.plot<-to.plot[stage1==stage2]
to.plot$group<-factor(to.plot$group,levels = c("intraembryo","interembryo"))

pdf("~/Desktop/NGS_postdoc/data/lineage_tracing/20190827.met.5kb.embryo.compare.box.pdf",width=5,height=3)
ggplot(to.plot[stage1==stage2],aes(stage1,cor,color=group))+
  geom_boxplot( outlier.size = 0.2,position = position_dodge(preserve = "single",width=0.9))+
  #geom_jitter()+
  ylim(0.1,0.7)+
  scale_color_manual(values=c("#DF8F44CC","#00A1D5CC"))+
  theme_bw()+
  labs(x="",y="Pairewised Spearman correaltions",title = "Methylome 5kb bin")+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=45,hjust=1),
    plot.title = element_text(hjust=0.5)
  )#+coord_fixed(ratio=4/0.6)
dev.off()

to.plot[,length(unique(cell1))]

to.plot[group=="intraembryo", mother_type:=ifelse( m_id1==m_id2,"same mother","different mother")]
to.plot$mother_type<-factor(to.plot$mother_type,levels = c("same mother","different mother"))
pdf("~/Desktop/NGS_postdoc/data/lineage_tracing/20190827.met.5kb.mother.compare.box.pdf",width=5,height=3)
ggplot(to.plot[group=="intraembryo" &stage1%in% c("2cell","4cell","L4cell","8cell")],aes(stage1,cor,color=mother_type))+
  geom_boxplot( outlier.size = 0.2,position = position_dodge(preserve = "single",width=0.9))+
  #geom_jitter()+
  ylim(0.1,0.7)+
  theme_bw()+
  scale_color_manual(values=c("#DC0000CC","#008280CC"))+
  labs(x="",y="Pairewised Spearman correaltions",title = "Methylome 5kb bin")+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=45,hjust=1),
    plot.title = element_text(hjust=0.5)
  )#+coord_fixed(ratio=4/0.6)
dev.off()

to.plot[group=="intraembryo", grand_type:=ifelse( g_id1==g_id2,"same grand","different grand")]
to.plot$grand_type<-factor(to.plot$grand_type,levels = c("same grand","different grand"))
pdf("~/Desktop/NGS_postdoc/data/lineage_tracing/20190827.met.5kb.grand.compare.box.pdf",width=5,height=3)
ggplot(to.plot[group=="intraembryo" &stage1%in% c("4cell","L4cell","8cell")],aes(stage1,cor,color=grand_type))+
  geom_boxplot( outlier.size = 0.2,position = position_dodge(preserve = "single",width=0.9))+
  #geom_jitter()+
  ylim(0.1,0.7)+
  scale_color_startrek()+
  theme_bw()+
  labs(x="",y="Pairewised Spearman correaltions",title = "Methylome 5kb bin")+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=45,hjust=1),
    plot.title = element_text(hjust=0.5)
  )#+coord_fixed(ratio=4/0.6)
dev.off()
### see 20200522.tracing.lineage.met.cortest.spearman.R  with 200 kb, cutoff 20 weight
to.plot<-rbind(to.plot[,.(Type=group,Stage=stage1,cor,Data="Met")],
      to.plot[group=="intraembryo" &stage1%in% c("2cell","4cell","L4cell","8cell"),.(Type=mother_type,Stage=stage1,cor,Data="Met")],
      to.plot[group=="intraembryo" &stage1%in% c("4cell","L4cell","8cell"),.(Type=grand_type,Stage=stage1,cor,Data="Met")]
)

####    acc NDR
## 20190717.clustering.R
correlation.test<-read.table("~/Desktop/NGS_postdoc/data/clustering/20190716.acc.ndr.spearman.correlations.tsv")
#rna.cor.dcast<-read.table("~/Desktop/NGS_postdoc/data/clustering/20190716.expr.tpm.spearman.correlations.tsv")

to.plot<-flattenCorrMatrix(correlation.test)%>% as.data.table %>% setnames(c("sample1","sample2","cor")) 
to.plot[,length(unique(sample1))]
#### 去除qc.acc.extra failed细胞

sample.anno<-fread("~/Desktop/NGS_postdoc/data/DNA.data/20190826.updated.sample.metadata.tsv")
sample.anno[qc.rna==T &qc.dna==T &qc.acc.extra==T,.N,stage]


to.plot<-merge(to.plot,sample.anno[qc.rna==T &qc.dna==T &qc.acc.extra==T,.(sample1=DNA.info,cell1=cellid,stage1=stage,embryo_id1=embryo_id,m_id1=mother_id,g_id1=g_id)],by="sample1")
to.plot<-merge(to.plot,sample.anno[qc.rna==T &qc.dna==T &qc.acc.extra==T,.(sample2=DNA.info,cell2=cellid,stage2=stage,embryo_id2=embryo_id,m_id2=mother_id,g_id2=g_id)],by="sample2")
sample.anno[cellid%in% to.plot$cell2,.N,sex]

unique(to.plot$embryo_id1)
to.plot[embryo_id1==embryo_id2 & embryo_id2=="X4cell2"]


unique(to.plot$stage1)
to.plot$stage1<-factor(to.plot$stage1,levels=c("zygote","2cell","4cell","L4cell","8cell","16cell", "ICM","TE")) #"ICM","TE"

to.plot[stage1==stage2, group:=ifelse(embryo_id1==embryo_id2,"intraembryo","interembryo")]
to.plot<-to.plot[stage1==stage2]
to.plot$group<-factor(to.plot$group,levels = c("intraembryo","interembryo"))

pdf("~/Desktop/NGS_postdoc/data/lineage_tracing/20190827.acc.ndr.embryo.compare.box.pdf",width=5,height=3)
ggplot(to.plot[stage1==stage2],aes(stage1,cor,color=group))+
  geom_boxplot( outlier.size = 0.2,position = position_dodge(preserve = "single",width=0.9))+
  #geom_jitter()+
  ylim(0,0.3)+
  scale_color_manual(values=c("#DF8F44CC","#00A1D5CC"))+
  theme_bw()+
  labs(x="",y="Pairewised Spearman correaltions",title = "Chromatin accessibility all NDR")+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=45,hjust=1),
    plot.title = element_text(hjust=0.5)
  )#+coord_fixed(ratio=4/0.6)
dev.off()

to.plot[,length(unique(cell1))]

to.plot[group=="intraembryo", mother_type:=ifelse( m_id1==m_id2,"same mother","different mother")]
to.plot$mother_type<-factor(to.plot$mother_type,levels = c("same mother","different mother"))
pdf("~/Desktop/NGS_postdoc/data/lineage_tracing/20190827.acc.ndr.mother.compare.box.pdf",width=5,height=3)
ggplot(to.plot[group=="intraembryo" &stage1%in% c("2cell","4cell","L4cell","8cell")],aes(stage1,cor,color=mother_type))+
  geom_boxplot( outlier.size = 0.2,position = position_dodge(preserve = "single",width=0.9))+
  #geom_jitter()+
  ylim(0,0.3)+
  theme_bw()+
  scale_color_manual(values=c("#DC0000CC","#008280CC"))+
  labs(x="",y="Pairewised Spearman correaltions",title = "Chromatin accessibility all NDR")+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=45,hjust=1),
    plot.title = element_text(hjust=0.5)
  )#+coord_fixed(ratio=4/0.6)
dev.off()

to.plot[group=="intraembryo", grand_type:=ifelse( g_id1==g_id2,"same grand","different grand")]
to.plot$grand_type<-factor(to.plot$grand_type,levels = c("same grand","different grand"))
pdf("~/Desktop/NGS_postdoc/data/lineage_tracing/20190827.acc.ndr.grand.compare.box.pdf",width=5,height=3)
ggplot(to.plot[group=="intraembryo" &stage1%in% c("4cell","L4cell","8cell")],aes(stage1,cor,color=grand_type))+
  geom_boxplot( outlier.size = 0.2,position = position_dodge(preserve = "single",width=0.9))+
  #geom_jitter()+
  ylim(0,0.3)+
  scale_color_startrek()+
  theme_bw()+
  labs(x="",y="Pairewised Spearman correaltions",title = "Chromatin accessibility all NDR")+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle=45,hjust=1),
    plot.title = element_text(hjust=0.5)
  )#+coord_fixed(ratio=4/0.6)
dev.off()
to.plot<-rbind(to.plot[,.(Type=group,Stage=stage1,cor,Data="Acc")],
               to.plot[group=="intraembryo" &stage1%in% c("2cell","4cell","L4cell","8cell"),.(Type=mother_type,Stage=stage1,cor,Data="Acc")],
               to.plot[group=="intraembryo" &stage1%in% c("4cell","L4cell","8cell"),.(Type=grand_type,Stage=stage1,cor,Data="Acc")]
)


flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}
