############################################################## 
##############################################################
## plot single genes
##############################################################

###### some NDR show overlap between distal and promoter, to assign the promoter NDR was based on the homer anno result, which use the center of NDR to the TSS, so it is possible some overlapped NDR assign to promoter NDR and distal NDR respectivly
# to plot single gene_NDR pairs, we need to remove the overlapped NDRs, to do so, we merge distal and promoter NDR, and then keep the max r pair NDRs
accrna.cor<-rbind(accrna.cor.p[,.(id,symbol,r,padj_fdr,type.ndr="promoter",sig.type,gene.type)],accrna.cor.d[,.(id,symbol,r,padj_fdr,type.ndr="distal",sig.type,gene.type)])
accrna.cor<-merge(ndr.anno[,.(id,chr,start,end,cpg_density)],accrna.cor,by="id") %>% setkey(chr,start,end) %>% .[,.(chr,start,end,id)]
fwrite(accrna.cor,"~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/plot.single.genes/20190813.all.ndrs.accrna.tsv",sep="\t",col.names = F)
#bedtools merge -i 20190813.all.ndrs.accrna.tsv -d 0  -c 4,4,4 -o distinct,collapse,count_distinct > 20190813.all.ndrs.accrna.merged.tsv
accrna.cor.merge<-fread("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/plot.single.genes/20190813.all.ndrs.accrna.merged.tsv",select = c(1:3)) %>% 
  setnames(c("chr","start","end")) %>% setkey(chr,start,end)
accrna.cor.merge[,nid:=paste0("nNDR",.I)]

accrna.cor<-rbind(accrna.cor.p[,.(id,symbol,r,padj_fdr,type.ndr="promoter",sig.type,gene.type)],accrna.cor.d[,.(id,symbol,r,padj_fdr,type.ndr="distal",sig.type,gene.type)])
accrna.cor<-merge(ndr.anno[,.(id,chr,start,end,cpg_density)],accrna.cor,by="id") %>% setkey(chr,start,end) #%>% .[,.(chr,start,end,id)]

test<-foverlaps(accrna.cor,accrna.cor.merge,nomatch = 0L)
dup.id<-test[,.N,nid]%>% .[N>1]
test[nid %in% dup.id$nid]
accrna.cor<-foverlaps(accrna.cor,accrna.cor.merge,nomatch = 0L)
accrna.cor<-accrna.cor[,.SD[which.max(r)],by="nid"] %>% .[,.(chr,start=i.start,end=i.end,id,symbol,r,padj_fdr,type.ndr,sig.type,gene.type)] # 114828 ->114782
fwrite(accrna.cor,"~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/plot.single.genes/20190813.all.ndrs.accrna.uniqued.in.both.distal.and.promoter.tsv",sep="\t")

##
test<-merge(ndr.anno[,.(chr,start,end,anno,id,type)],accrna.cor.p[sig.type=="pos.sig",.(id,symbol=paste0(symbol,"_",round(r*100)/100))],by="id") %>% setkey(chr,start,end) 
test<-merge(ndr.anno[,.(chr,start,end,anno,id,type)],accrna.cor.d[sig.type=="pos.sig",.(id,symbol=paste0(symbol,"_",round(r*100)/100))],by="id") %>% setkey(chr,start,end) 
test[symbol=="Pou5f1_0.48"]
gene.anno[symbol=="Pou5f1"]

sample.anno<-fread("~/Desktop/NGS_postdoc/data/DNA.data/20190524.updated.cell.metadata.tsv")
sample.anno$stage<-factor(sample.anno$stage,levels =c("zygote","2cell","4cell","L4cell","8cell","16cell","ICM","TE"))

acc.sc<-fread("~/Desktop/NGS_postdoc/data/DNA.data/201903.parsed.epi/20190529.sc.acc.NDRs.hiAccRNA.w.N.Nacc.tsv.gz")
acc.sc<-merge(acc.sc,sample.anno[qc.rna==T &qc.dna==T & qc.acc.extra==T,.(sample=DNA.info,cellid,stage,DC.order)],by="sample")
acc.sc<-acc.sc[weight>=3]

accrna.cor<-fread("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/plot.single.genes/20190813.all.ndrs.accrna.uniqued.in.both.distal.and.promoter.tsv")

fwrite(accrna.cor[sig.type=="pos.sig",.(chr,start,end,id=paste0(symbol,"_",round(r*100)/100))],"~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/plot.single.genes/20190813.ndrs.pos.sig.accrna.bed",sep="\t",col.names = F)
fwrite(accrna.cor[,.(chr,start,end,id=paste0(symbol,"_",round(r*100)/100))],"~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/plot.single.genes/20190813.ndrs.all.uniqued.accrna.bed",sep="\t",col.names = F)


rna.tpm<-fread("~/Desktop/NGS_postdoc/data/RNA.expr.mat/20190519.tech.rep.collapsed.RNA.matrix.tpm.tsv")
gene.anno<-fread("~/Desktop/NGS_postdoc/data/DNA.data/ref/wy.gencode.gene.metadata.tsv")%>% .[,.(gene=ens_ID,symbol)]
#sample.anno<-fread("~/Desktop/NGS_postdoc/data/DNA.data/20190524.updated.cell.metadata.tsv") %>% .[qc.rna==T]
unique(sample.anno$stage)
rna.tpm<-merge(rna.tpm,sample.anno[qc.rna==T &qc.dna==T &qc.acc.extra==T,.(cellid,stage)],by="cellid")
rna.tpm[,expr:=log2(tpm+1)]
rna.tpm<-merge(rna.tpm,gene.anno,by="gene")


accrna.cor[gene.type=="TE.high" & sig.type=="pos.sig" &r>0.5,unique(symbol)]

 for (the.gene in c( "Cdx2"  , "Nanog","Pou5f1") ) {
#accrna.cor[symbol==the.gene]
to.plot.sg<-merge(accrna.cor[symbol==the.gene,.(id,symbol,r,fdr=-log10(padj_fdr),sig.type,type.ndr)],acc.sc[,.(id,mean_rate,weight,cellid,stage,DC.order)],by="id")
to.plot.sg<-merge(to.plot.sg,ndr.anno[,.(chr,start,end,id)],by="id") %>% setkey(chr,start,end)
to.plot.sg[,plot.id:=paste0(chr,":",start,"-",end)]
to.plot.sg[sig.type=="pos.sig",unique(plot.id)]
#fwrite(to.plot.sg[,.(chr,start,end,symbol=paste0(symbol,"_",round(r*100)/100))]%>% unique(), sprintf("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/plot.single.genes/%s.all.bed",the.gene),sep="\t",col.names = F)
#fwrite(to.plot.sg[sig.type=="pos.sig",.(chr,start,end,symbol=paste0(symbol,"_",round(r*100)/100))]%>% unique(), sprintf("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/plot.single.genes/%s.pos.sig.bed",the.gene),sep="\t",col.names = F)

to.plot.sg.matrix<-dcast(to.plot.sg,cellid~plot.id,value.var = "mean_rate")%>% as.data.frame() %>% tibble::column_to_rownames("cellid")
row.anno<-sample.anno[cellid%in% rownames(to.plot.sg.matrix),.(cellid,stage,DC.order)]%>% setkey(stage,DC.order) %>% .[,DC.order:=NULL] %>% as.data.frame() %>% tibble::column_to_rownames("cellid")
#dev.off()
library(RColorBrewer)
pdf(sprintf("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/plot.single.genes/20190813.acc.%s.pdf",the.gene),width=8,height=5)
p1<-pheatmap::pheatmap(to.plot.sg.matrix[rownames(row.anno),], cluster_rows = F,cluster_cols = F,annotation_row = row.anno,show_rownames = F,
                   border_color = NA,na_col = "white",cellwidth = 5,cellheight = 1,
                   breaks = seq(0,100,by=5), fontsize =5,
                   color = colorRampPalette(c("#006699","#CCCCFF","#FFCC00","red"),bias=0.5)(length( seq(0,100,by=5))), #
                   #color = colorRampPalette((brewer.pal(n = 9, name ="YlGnBu")))(length( seq(0,100,by=5))), 
                   annotation_colors=list(stage = c( zygote="#7FC97F",`2cell`= "#BEAED4", `4cell`= "#FDC086",L4cell= "#808000",
                                                     `8cell`= "#386CB0",`16cell`= "#650136", ICM= "#E64B35CC",TE= "#4DBBD5CC"))
)
print(p1)
dev.off()
 }

to.plot.sg.matrix.r<-dcast(to.plot.sg[,.(plot.id,r,plots="1")] %>% unique,plots~plot.id,value.var = "r")%>% as.data.frame() %>% tibble::column_to_rownames("plots")
min(to.plot.sg.matrix.r)
max(to.plot.sg.matrix.r)

pdf(sprintf("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/plot.single.genes/20190813.acc.%s.accrna.cor.pdf",the.gene),width=8,height=3)
pheatmap::pheatmap(to.plot.sg.matrix.r, cluster_rows = F,cluster_cols = F,show_rownames = F,
                   border_color = "grey90",na_col = "grey90",cellwidth = 5,cellheight = 5,
                   breaks = seq(-0.7,0.7,by=0.1),fontsize =5,
                   color = colorRampPalette(c("royalblue","white","orangered"))(length( seq(-0.7,0.7,by=0.1)))
)
dev.off()

to.plot.sg.matrix.fdr<-dcast(to.plot.sg[,.(plot.id,fdr,plots="1")] %>% unique,plots~plot.id,value.var = "fdr")%>% as.data.frame() %>% tibble::column_to_rownames("plots")
to.plot.sg.matrix.fdr[to.plot.sg.matrix.fdr<=1]<-0
min(to.plot.sg.matrix.fdr)
max(to.plot.sg.matrix.fdr)
pdf(sprintf("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/plot.single.genes/20190813.acc.%s.accrna.cor.fdr.pdf",the.gene),width=8,height=3)
pheatmap::pheatmap(to.plot.sg.matrix.fdr, cluster_rows = F,cluster_cols = F,show_rownames = F,
                   border_color = "grey90",na_col = "grey90",cellwidth = 5,cellheight = 5,
                   breaks = seq(0,5,by=0.5),fontsize = 5,
                   color = colorRampPalette(c("white","green"))(length( seq(0,5,by=0.5)))
)
dev.off()



to.plot.sg.matrix.expr<-dcast(rna.tpm[symbol==the.gene,.(cellid,expr,plots="1")],cellid~plots,value.var = "expr") %>% as.data.frame() %>% tibble::column_to_rownames("cellid")
max(to.plot.sg.matrix.expr)
min(to.plot.sg.matrix.expr)
pdf(sprintf("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/plot.single.genes/20190813.acc.%s.expr.pdf",the.gene),width=8,height=5)
pheatmap::pheatmap(to.plot.sg.matrix.expr[rownames(row.anno),],scale="column", cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,
                   border_color = NA,na_col = "grey90",cellwidth = 10,cellheight = 1,
                   breaks = seq(-2,2,by=0.2), fontsize =5,
                   color = colorRampPalette(c("black","yellow"))(length( seq(-2,2,by=0.2))),
                   annotation_colors=list(stage = c( zygote="#7FC97F",`2cell`= "#BEAED4", `4cell`= "#FDC086",L4cell= "#808000",
                                                     `8cell`= "#386CB0",`16cell`= "#650136", ICM= "#E64B35CC",TE= "#4DBBD5CC"))
)
dev.off()


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## plot single NDR_gene cor   x=expr y=acc r p, chr:start-end
to.plot<-accrna.cor[symbol%in% c("Nanog","Pou5f1","Cdx2") &sig.type=="pos.sig"]
to.plot[,labels:=paste0(symbol,"_",chr,":",start,"-",end,"_",round(r*100)/100)]

to.plot<-merge(to.plot,acc.sc[,.(id,mean_rate,weight,cellid,stage)],by="id") %>% merge(rna.tpm[,.(symbol,cellid,expr)],by=c("symbol","cellid"))
i="Cdx2_chr5:148023421-148023620_0.29"

for( i in to.plot[symbol=="Pou5f1",unique(labels)]) {
  
  x.len=to.plot[labels==i,max(expr)-min(expr)]
  y.len=to.plot[labels==i,max(mean_rate)-min(mean_rate)]
  r.value=to.plot[labels==i,unique(r)]
  pdf(sprintf("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/plot.single.genes/20190819.scater.%s.pdf",i),width = 4,height=3)
  p1<-ggplot(to.plot[labels==i],aes(expr,mean_rate))+
    geom_point(size=1,alpha=0.8,aes(color=stage))+
    scale_color_manual(values = c(  zygote="#7FC97F",`2cell`= "#BEAED4", `4cell`= "#FDC086",L4cell= "#808000",
                                    `8cell`= "#386CB0",`16cell`= "#650136", ICM= "#E64B35CC",TE= "#4DBBD5CC"))+ 
    theme_bw()+
    scale_y_continuous(breaks = c(0,50,100))+
    # scale_x_continuous(breaks = c(0,5,10))+
    scale_x_continuous(breaks = c(8,10,12))+
    geom_smooth(method = lm,se = F)+
    #geom_abline(slope=r.value)+
    labs(title = i,x="RNA level",y="Acc level")+
    coord_fixed(ratio = x.len/y.len)
  print(p1)
  dev.off()
}