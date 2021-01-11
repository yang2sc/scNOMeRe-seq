

##### figure 5e
##### compare met and expr in ACC/RNA inferred CREs 
acc.met.rna.cor.d<-fread("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/20190810.distal.NDR_accmetrnaTPM.wtd.cor.final.tsv")
acc.met.rna.cor.p<-fread("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/20190810.proximal.NDR_accmetrnaTPM.wtd.cor.final.tsv")


tmp<-rbind(acc.met.rna.cor.d[,NDR.type:="distal"], acc.met.rna.cor.p[,NDR.type:="promoter"])

for (sig.types in c("pos.sig","neg.sig") ) {
  for (ndr.types in c("distal", "promoter") ) {
    to.plot<-tmp[acc.sig.type==sig.types & NDR.type==ndr.types]
    setorder(to.plot,-met.p)
    y.limit<-to.plot[,max(-log10(met.p))]
    all<-nrow(to.plot)
    all.gene<-to.plot[,length(unique(symbol))]
    
    all.icm<-nrow(to.plot[gene.type=="ICM.high"])
    all.te<-nrow(to.plot[gene.type=="TE.high"])
    all.gene.icm<-to.plot[gene.type=="ICM.high",length(unique(symbol))]
    all.gene.te<-to.plot[gene.type=="TE.high",length(unique(symbol))]
    
    negative_hits.icm<-to.plot[met.sig.type== "neg.sig" &gene.type=="ICM.high",symbol]
    negative.gene.icm<-to.plot[met.sig.type== "neg.sig"&gene.type=="ICM.high",length(unique(symbol))]
    positive_hits.icm<-to.plot[met.sig.type== "pos.sig"&gene.type=="ICM.high",symbol]
    positive.gene.icm<-to.plot[met.sig.type== "pos.sig"&gene.type=="ICM.high",length(unique(symbol))]
    
    negative_hits.te<-to.plot[met.sig.type== "neg.sig"&gene.type=="TE.high",symbol]
    negative.gene.te<-to.plot[met.sig.type== "neg.sig" &gene.type=="TE.high",length(unique(symbol))]
    positive_hits.te<-to.plot[met.sig.type== "pos.sig" &gene.type=="TE.high",symbol]
    positive.gene.te<-to.plot[met.sig.type== "pos.sig" &gene.type=="TE.high",length(unique(symbol))]
    
    
    pdf(paste0("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/accrna.metrna.defined.cres/20190915.icm.te.accrna.",sig.types,".metrna.",ndr.types,".wtd.cor.vocalno.pdf"),width=6,height=3)
    p1<-ggplot(to.plot,aes(x=met.r ,-log10(met.p)))+ # y=(-log10(acc.p))))+
      geom_segment(aes(x=0, xend=0, y=0, yend=y.limit*0.8), color="grey",linetype="dashed") +
      #  geom_vline(xintercept = 0,linetype="dashed")+
      geom_point(aes(color=met.sig.type),size=0.2,alpha=0.5)+
      scale_color_manual(values=c("#386CB0","grey","#E64B35CC")) +  #"#7FC97F","#386CB0","#BEAED4",
      labs(title=paste0(sig.types,"_Acc/RNA_",ndr.types))+
      xlab("Met/RNA correlation")+ #Weighted Pearson correlation
      ylab("-log10(p)")+
      xlim(-1,1)+
      #  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
      # scale_y_continuous(limits=c(0,15)) +
      annotate("text", x=0, y=y.limit*0.95, size=4, label=sprintf("All\n%d/%d\n%d/%d", all.icm,all.gene.icm,all.te,all.gene.te)) +
      annotate("text", x=-0.5, y=y.limit*0.95, size=4, label=sprintf("(-)\n%d/%d\n%d/%d",length(negative_hits.icm),negative.gene.icm,length(negative_hits.te),negative.gene.te)) +
      annotate("text", x=0.5, y=y.limit*0.95, size=4, label=sprintf("(+)\n%d/%d\n%d/%d",length(positive_hits.icm),positive.gene.icm,length(positive_hits.te),positive.gene.te)) +
      annotate("text", x= 0, y=y.limit*0.7, size=4, label=sprintf("\nICM.high (NDR/genes)\nTE.high (NDR/genes)")) +
      theme_classic()+
      guides(color=F)+
      theme(
        legend.position = "bottom",
        axis.title = element_text(size=16), # family="Arial"
        axis.text = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_blank()#,
        #  legend.direction = "vertical"
      )
    print(p1)
    dev.off()
    
  }}

