# figure 2 a-g
library(data.table)
library(purrr)
library(furrr)
library(dplyr)
library(plyr)

sample.anno<-fread("~/Desktop/NGS_postdoc/data/DNA.data/20190524.updated.cell.metadata.tsv")
##### add lineage tracing information
sample.anno$id<-do.call(rbind,strsplit(sample.anno$cellid,"[_]"))[,2]

sample.anno[stage %in% c("2cell"),mother_id:=embryo_id]
sample.anno[stage %in% c("4cell","L4cell","8cell") & id %in% c(1,2),mother_id:=paste0(embryo_id,"_a")]
sample.anno[stage %in% c("4cell","L4cell","8cell") & id %in% c(3,4),mother_id:=paste0(embryo_id,"_b")]
sample.anno[stage %in% c("8cell") & id %in% c(5,6),mother_id:=paste0(embryo_id,"_c")]
sample.anno[stage %in% c("8cell") & id %in% c(7,8),mother_id:=paste0(embryo_id,"_d")]

sample.anno[stage %in% c("4cell","L4cell") ,g_id:=embryo_id]
sample.anno[embryo_id =="X8cell1" & id %in% c(1:4),g_id:=paste0(embryo_id,"_e")]
sample.anno[embryo_id =="X8cell1" & id %in% c(5:8),g_id:=paste0(embryo_id,"_f")]

sample.anno[embryo_id =="X8cell2" & id %in% c(1,2,7,8),g_id:=paste0(embryo_id,"_e")]
sample.anno[embryo_id =="X8cell2" & id %in% c(3:6),g_id:=paste0(embryo_id,"_f")]

sample.anno[embryo_id =="X8cell3" & id %in% c(1,2,5,6),g_id:=paste0(embryo_id,"_e")]
sample.anno[embryo_id =="X8cell3" & id %in% c(3,4,7,8),g_id:=paste0(embryo_id,"_f")]

sample.anno[embryo_id =="X8cell4" & id %in% c(1,2,5,6),g_id:=paste0(embryo_id,"_e")]
sample.anno[embryo_id =="X8cell4" & id %in% c(3,4,7,8),g_id:=paste0(embryo_id,"_f")]
sample.anno$id<-NULL

sample.anno[!is.na(mother_id)]
fwrite(sample.anno,'~/Desktop/NGS_postdoc/data/DNA.data/20190826.updated.sample.metadata.tsv',sep="\t")


met.1mb<-fread("~/Desktop/NGS_postdoc/data/DNA.data/201903.parsed.epi/20190423.parsed.met.1M.bin.tsv") # parsed in 20190128.cnv.R
bin.metadata<-fread("~/Desktop/NGS_postdoc/data/DNA.data/201903.parsed.epi/20190423.parsed.met.1M.bin.parse.metadata.tsv")
sample.anno<-fread("~/Desktop/NGS_postdoc/data/DNA.data/20190826.updated.sample.metadata.tsv")
met.1mb<-merge(met.1mb,sample.anno[qc.rna==T &qc.dna==T & !sex=="GG",.(sample=DNA.info,cellid,sex,embryo_id,mother_id,g_id,stage)],by="sample")
summary(met.1mb$weight)
met.1mb[weight<3]

# cal zscore by the embryo_id and bin id
tmp<-met.1mb[,sum(weight),cellid] %>% setorder(V1)
head(tmp,20)
met.1mb[weight>=3,zscore:=(mean_rate-mean(mean_rate))/sd(mean_rate),.(embryo_id,id)]

met.1mb.dcast<-dcast(met.1mb[!is.na(zscore) & !stage=="zygote"],id~cellid, value.var = "zscore")%>% as.data.frame() %>% tibble::column_to_rownames("id")
#met.1mb_correlations<-cor(met.1mb.dcast,use="pairwise.complete.obs",method="spearman")
met.1mb_correlations<-cor(met.1mb.dcast,use="pairwise.complete.obs",method="pearson")

test<-melt(met.1mb_correlations)%>% as.data.table %>% setnames(c("cell1","cell2","cor")) 

test<-merge(test,sample.anno[,.(cell1=cellid,stage1=stage,embryo_id1=embryo_id,m_id1=mother_id,g_id1=g_id)],by="cell1", all.x=T)
test<-merge(test,sample.anno[,.(cell2=cellid,stage2=stage,embryo_id2=embryo_id,m_id2=mother_id,g_id2=g_id)],by="cell2", all.x=T)
unique(test$embryo_id1)
test[embryo_id1==embryo_id2 & embryo_id2=="X4cell2"]

to.plot<-dcast(test[embryo_id1==embryo_id2 & embryo_id1=="X4cell9"],cell1~cell2,value.var = "cor") %>% as.data.frame() %>% tibble::column_to_rownames("cell1")

pdf("~/Desktop/NGS_postdoc/data/lineage_tracing/20190827.tracing.1M.bin.pearson.cor.4cell9.heatmap.pdf")
pheatmap::pheatmap(to.plot,breaks = seq(-1,1,by=0.1),# clustering_method = "ward.D2",
                   cellwidth = 20,cellheight = 20,
                   display_numbers = T,number_format = "%.2f",number_color = "white",
                   angle_col = "45", fontsize = 10,
                   color = colorRampPalette(c("darkblue","grey","orangered"))(length( seq(-1,1,by=0.1))) )
dev.off()

to.plot<-dcast(test[embryo_id1==embryo_id2 & embryo_id1=="L4cell6"],cell1~cell2,value.var = "cor") %>% as.data.frame() %>% tibble::column_to_rownames("cell1")

pdf("~/Desktop/NGS_postdoc/data/lineage_tracing/20190827.tracing.1M.bin.pearson.cor.L4cell6.heatmap.pdf")
pheatmap::pheatmap(to.plot,breaks = seq(-1,1,by=0.1),# clustering_method = "ward.D2",
                   cellwidth = 20,cellheight = 20,
                   display_numbers = T,number_format = "%.2f",number_color = "white",
                   angle_col = "45", fontsize = 10,
                   color = colorRampPalette(c("darkblue","grey","orangered"))(length( seq(-1,1,by=0.1))) )
dev.off()

to.plot<-dcast(test[embryo_id1==embryo_id2 & embryo_id1=="X8cell4"],cell1~cell2,value.var = "cor") %>% as.data.frame() %>% tibble::column_to_rownames("cell1")

pdf("~/Desktop/NGS_postdoc/data/lineage_tracing/20190827.tracing.1M.bin.pearson.cor.8cell4.heatmap.pdf")
pheatmap::pheatmap(to.plot,breaks = seq(-1,1,by=0.1),# clustering_method = "ward.D2",
                   cellwidth = 20,cellheight = 20,
                   display_numbers = T,number_format = "%.2f",number_color = "white",
                   angle_col = "45", fontsize = 10,
                   color = colorRampPalette(c("darkblue","grey","orangered"))(length( seq(-1,1,by=0.1))) )
dev.off()

#density plot

flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}
to.plot<-flattenCorrMatrix(met.1mb_correlations)%>% as.data.table %>% setnames(c("cell1","cell2","cor")) 
to.plot<-merge(to.plot,sample.anno[,.(cell1=cellid,stage1=stage,embryo_id1=embryo_id,m_id1=mother_id,g_id1=g_id)],by="cell1", all.x=T)
to.plot<-merge(to.plot,sample.anno[,.(cell2=cellid,stage2=stage,embryo_id2=embryo_id,m_id2=mother_id,g_id2=g_id)],by="cell2", all.x=T)
unique(to.plot$embryo_id1)
to.plot[embryo_id1==embryo_id2 & embryo_id2=="X4cell2"]

to.plot$stage1<-factor(to.plot$stage1,levels=c("zygote","2cell","4cell","L4cell","8cell","16cell","ICM","TE"))
to.plot<-to.plot[embryo_id1==embryo_id2 &stage1 %in% c("4cell","L4cell","8cell")]#,"16cell","ICM","TE")]


pdf("~/Desktop/NGS_postdoc/data/lineage_tracing/20190827.zscore.1mbin.met.pearson.density.pdf",width = 4,height=6)
ggplot(to.plot[embryo_id1==embryo_id2 ], aes(x=cor, fill=stage1)) +
  geom_density(alpha=0.8,color=NA)+
  facet_grid(stage1~.)+
  scale_fill_jco()+
  scale_y_continuous(breaks = c(0,1,2))+
  xlab("Pairewise correlation of blastomeres")+
  ylab("Density")+
  theme_bw()+
  theme(
    #    panel.grid.minor.y = element_blank(),
    #    panel.grid.major.y = element_blank(),
    strip.background = element_blank()
  )
dev.off()

to.plot<-to.plot[embryo_id1==embryo_id2 &stage1 %in% c("4cell","L4cell","8cell")]#,"16cell","ICM","TE")]

pdf("~/Desktop/NGS_postdoc/data/lineage_tracing/20190827.zscore.1mbin.met.pearson.density.4.8cell.pdf",width = 4,height=4)
ggplot(to.plot[embryo_id1==embryo_id2 ], aes(x=cor, fill=stage1)) +
  geom_density(alpha=0.8,color=NA)+
  facet_grid(stage1~.)+
  #scale_fill_jco()+
  scale_y_continuous(breaks = c(0,1,2))+
  xlab("Pairewise correlation of blastomeres")+
  ylab("Density")+
  theme_bw()+
  theme(
    #    panel.grid.minor.y = element_blank(),
    #    panel.grid.major.y = element_blank(),
    strip.background = element_blank()
  )
dev.off()

#############   validation 
########################################## FITC and validation BS-seq
bs.seq.sample.info<-fread("~/Desktop/NGS_postdoc/data/lineage_tracing/20190424.bs.seq.sample.info.txt")
bs.seq.sample.info$embryo_id<-do.call(rbind,strsplit(bs.seq.sample.info$cellid,"[_]"))[,2]
bs.seq.sample.info$type<-do.call(rbind,strsplit(bs.seq.sample.info$cellid,"[_]"))[,1]
bs.seq.sample.info[,stage:=substring(embryo_id,1,5)]
bs.seq.sample.info[,embryo_id:=paste0(type,"_",embryo_id)]

fwrite(bs.seq.sample.info,"~/Desktop/NGS_postdoc/data/DNA.data/20190424.bsseq.sample.info.tsv",sep="\t")
bs.seq.sample.info<-fread("~/Desktop/NGS_postdoc/data/DNA.data/20190424.bsseq.sample.info.tsv")
bs.seq.1m.bin<-fread("~/Desktop/NGS_postdoc/data/DNA.data/201903.parsed.epi/20190423.parsed.met.1M.bin.bs_seq.tsv")
bs.seq.1m.bin<-merge(bs.seq.1m.bin,bs.seq.sample.info,by.x="sample",by.y="DNA.info")
summary(bs.seq.1m.bin$weight)
bs.seq.1m.bin[weight<10]
bs.seq.1m.bin[,sum(weight),cellid]#FITC_8cell2_7    8068
bs.seq.1m.bin<-bs.seq.1m.bin[!cellid=="FITC_8cell2_7"]

bs.seq.1m.bin[,zscore:=(mean_rate-mean(mean_rate))/sd(mean_rate),.(embryo_id,id)]

bs.seq.1m.bin<-bs.seq.1m.bin[!is.na(zscore)]
bs.seq.1m.bin.dcast<-dcast(bs.seq.1m.bin,id~cellid, value.var = "zscore")%>% as.data.frame() %>% tibble::column_to_rownames("id")
#met.1mb_correlations<-cor(met.1mb.dcast,use="pairwise.complete.obs",method="spearman")
bs.seq_correlations<-cor(bs.seq.1m.bin.dcast,use="pairwise.complete.obs",method="pearson")

test<-flattenCorrMatrix(bs.seq_correlations)%>% as.data.table %>% setnames(c("cell1","cell2","cor")) 
test<-melt(bs.seq_correlations)%>% as.data.table %>% setnames(c("cell1","cell2","cor")) 

test<-merge(test,bs.seq.sample.info[,.(cell1=cellid,stage1=stage,embryo_id1=embryo_id)],by="cell1", all.x=T)
test<-merge(test,bs.seq.sample.info[,.(cell2=cellid,stage2=stage,embryo_id2=embryo_id)],by="cell2", all.x=T)
unique(test$embryo_id1)

pdf("~/Desktop/NGS_postdoc/data/lineage_tracing/20190827.zscore.1mbin.met.bs.seq.8cell1.pearson.heatmap.pdf")
to.plot<-dcast(test[embryo_id1==embryo_id2 & embryo_id1=="FITC_8cell1" ],cell1~cell2,value.var = "cor") %>% as.data.frame() %>% tibble::column_to_rownames("cell1")
pheatmap::pheatmap(to.plot,breaks = seq(-1,1,by=0.1),# clustering_method = "ward.D2",
                   cellwidth = 20,cellheight = 20,border_color = "white",
                   display_numbers = T,number_format = "%.2f",number_color = "white",
                   angle_col = "45", fontsize = 10,
                   color = colorRampPalette(c("darkblue","grey","orangered"))(length( seq(-1,1,by=0.1))) )
dev.off()

pdf("~/Desktop/NGS_postdoc/data/lineage_tracing/20190827.zscore.1mbin.met.bs.seq.8cell2.pearson.heatmap.pdf")
to.plot<-dcast(test[embryo_id1==embryo_id2 & embryo_id1=="FITC_8cell2" ],cell1~cell2,value.var = "cor") %>% as.data.frame() %>% tibble::column_to_rownames("cell1")
pheatmap::pheatmap(to.plot,breaks = seq(-1,1,by=0.1),# clustering_method = "ward.D2",
                   cellwidth = 20,cellheight = 20,border_color = "white",
                   display_numbers = T,number_format = "%.2f",number_color = "white",
                   angle_col = "45", fontsize = 10,
                   color = colorRampPalette(c("darkblue","grey","orangered"))(length( seq(-1,1,by=0.1))) )
dev.off()

pdf("~/Desktop/NGS_postdoc/data/lineage_tracing/20190827.zscore.1mbin.met.bs.seq.4cell1.pearson.heatmap.pdf")
to.plot<-dcast(test[embryo_id1==embryo_id2 & embryo_id1=="FITC_4cell1" ],cell1~cell2,value.var = "cor") %>% as.data.frame() %>% tibble::column_to_rownames("cell1")
pheatmap::pheatmap(to.plot,breaks = seq(-1,1,by=0.1),# clustering_method = "ward.D2",
                   cellwidth = 20,cellheight = 20,border_color = "white",
                   display_numbers = T,number_format = "%.2f",number_color = "white",
                   angle_col = "45", fontsize = 10,
                   color = colorRampPalette(c("darkblue","grey","orangered"))(length( seq(-1,1,by=0.1))) )
dev.off()

