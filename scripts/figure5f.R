# figure 5f
library(data.table)
library(purrr)
library(furrr)
library(dplyr)
library(plyr)

## compare met and expr in ACC/RNA inferred CREs 
acc.met.rna.cor.d<-fread("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/20190810.distal.NDR_accmetrnaTPM.wtd.cor.final.tsv")
acc.met.rna.cor.p<-fread("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/20190810.proximal.NDR_accmetrnaTPM.wtd.cor.final.tsv")

accmetrna.cor<-rbind(acc.met.rna.cor.d[,NDR.type:="distal"], acc.met.rna.cor.p[,NDR.type:="promoter"])

diff.acc<-fread("~/Desktop/NGS_postdoc/data/te.icm/diff.acc.met/20190625.ICM_vs_TE.acc.NDR.diff.test.tsv") %>% .[!is.na(diff)] #20190625.ICM_vs_TE.acc.NDR.value.tsv
diff.met<-fread("~/Desktop/NGS_postdoc/data/te.icm/diff.acc.met/20190625.ICM_vs_TE.met.NDR.diff.test.tsv") %>% .[!is.na(diff)] #20190625.ICM_vs_TE.acc.NDR.value.tsv

accmetrna.diff<-merge(accmetrna.cor[!acc.sig.type=="not.sig",.(id,symbol,NDR.type,sig.type,gene.type,acc.sig.type,met.sig.type)],diff.acc[,.(id,diff.acc=diff,icm.acc=prop1,te.acc=prop2,acc.sig=sig)],by="id") %>% merge(diff.met[,.(id,diff.met=diff,icm.met=prop1,te.met=prop2,met.sig=sig)],by="id")

accmetrna.diff$gene.type<-factor(accmetrna.diff$gene.type,levels = c("ICM.high","TE.high"))
accmetrna.diff$acc.sig.type<-factor(accmetrna.diff$acc.sig.type,levels = c("pos.sig","neg.sig"))
accmetrna.diff$met.sig.type<-factor(accmetrna.diff$met.sig.type,levels = c("pos.sig","not.sig","neg.sig"))
accmetrna.diff$NDR.type<-factor(accmetrna.diff$NDR.type,levels = c("promoter","distal"))

accmetrna.diff[sig.type=="ACC.pos.sig_MET.neg.sig" &gene.type=="ICM.high" &met.sig==TRUE &icm.met<te.met,unique(symbol)]
accmetrna.diff[sig.type=="ACC.pos.sig_MET.neg.sig" &gene.type=="ICM.high" &met.sig==TRUE &icm.met<te.met]
accmetrna.diff[sig.type=="ACC.pos.sig_MET.neg.sig" &gene.type=="TE.high" &met.sig==TRUE &icm.met>te.met,unique(symbol)]

#accrna.diff[,plot.type:=paste0(gene.type,"_Acc.",acc.sig,"_Met.",met.sig)]
pdf("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/accrna.metrna.defined.cres/20190915.accmetrna.cor.diff.acc.diff.met.distal.pdf",width = 6,height = 4)
ggplot(accmetrna.diff[NDR.type=="distal" ],aes(diff.acc,diff.met,color=gene.type))+
  geom_point(alpha=0.5,size=0.2)+
scale_color_manual(values=c(ICM.high= "#E64B35CC",TE.high= "#4DBBD5CC"))+
  facet_grid( met.sig.type~acc.sig.type)+
  geom_vline(xintercept = 0,linetype="dashed",size=0.5,alpha=0.5,color="black")+
  geom_hline(yintercept = 0,linetype="dashed",size=0.5,alpha=0.5,color="black")+
  labs(x="Differential Acc level\n(TE-ICM)",y="Differential Met level\n(TE-ICM)", title = "Distal CREs")+
  scale_y_continuous(breaks = c(0.5, 0, -0.5),limits =  c(-0.5,0.5))+
  xlim(-0.8, 0.8)+
  theme_bw()+
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    axis.title = element_text(size=12),
    axis.text = element_text(size=10),
    strip.text = element_text(size=12),
    plot.title = element_text(hjust = 0.5)
  )+ coord_fixed(ratio=1)
dev.off()

pdf("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/accrna.metrna.defined.cres/20190915.accmetrna.cor.diff.acc.diff.met.promoter.pdf",width = 6,height = 4)
ggplot(accmetrna.diff[NDR.type=="promoter" ],aes(diff.acc,diff.met,color=gene.type))+
  geom_point(alpha=0.5,size=0.2)+
  scale_color_manual(values=c(ICM.high= "#E64B35CC",TE.high= "#4DBBD5CC"))+
  facet_grid( met.sig.type~acc.sig.type)+
  geom_vline(xintercept = 0,linetype="dashed",size=0.5,alpha=0.5,color="black")+
  geom_hline(yintercept = 0,linetype="dashed",size=0.5,alpha=0.5,color="black")+
  labs(x="Differential Acc level\n(TE-ICM)",y="Differential Met level\n(TE-ICM)", title = "Promoter CREs")+
  scale_y_continuous(breaks = c(0.5, 0, -0.5),limits =  c(-0.5,0.5))+
  xlim(-0.8, 0.8)+
  theme_bw()+
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    axis.title = element_text(size=12),
    axis.text = element_text(size=10),
    strip.text = element_text(size=12),
    plot.title = element_text(hjust = 0.5)
  )+ coord_fixed(ratio=1)
dev.off()
