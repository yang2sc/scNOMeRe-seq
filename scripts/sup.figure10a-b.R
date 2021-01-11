#extended data 9a-d
library(ggpubr)
#### plot acc met
#### met

tmp<-list.files("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/random.set/met.level",pattern=".NDR3000.w.random.tsv",full.names = T)
tmp<-lapply(tmp, function(f) fread(f) %>% .[,stage:= basename(f) %>% 
                                              sub("20190813.local.met.level.","",.) %>% 
                                              sub(".NDR3000.w.random.tsv","",.)]) %>% rbindlist()

#### acc 
tmp<-list.files("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/random.set/acc.level",pattern=".NDR3000.w.random.tsv",full.names = T)
tmp<-lapply(tmp, function(f) fread(f) %>% .[,stage:= basename(f) %>% 
                                              sub("20190813.local.acc.level.","",.) %>% 
                                              sub(".NDR3000.w.random.tsv","",.)]) %>% rbindlist()


tmp$stage<-factor(tmp$stage,levels = c("zygote","2cell","4cell","L4cell","8cell","16cell","ICM","TE"))
tmp[,sig.type:=gene.type]
tmp$sig.type<-do.call(rbind,strsplit(tmp$sig.type,"_"))[,2]
tmp$anno<-do.call(rbind,strsplit(tmp$gene.type,"_"))[,3]
tmp$gene.type<-do.call(rbind,strsplit(tmp$gene.type,"_"))[,1]
tmp$sig.type<-factor(tmp$sig.type,levels = c('pos.sig',"neg.sig"))
tmp$anno<-factor(tmp$anno,levels=c("CRE","random"))
tmp$gene.type<-factor(tmp$gene.type,levels = c("ICM.high","TE.high"))
tmp$ndr.types<-factor(tmp$ndr.types,levels=c("promoter","distal"))

#tmp$sig.type<-factor(tmp$sig.type,levels = c( "ICM.high_CRE","TE.high_CRE","ICM.high_random","TE.high_random"))
#tmp$ndr.type<-factor(tmp$ndr.type,levels = c("promoter","distal"))
#tmp$gene.type<-factor(tmp$gene.type,levels = c( "ICM.high_CRE","TE.high_CRE","ICM.high_random","TE.high_random"))


# met
sample.anno<-fread("~/Desktop/NGS_postdoc/data/DNA.data/20190524.updated.cell.metadata.tsv")
tmp<-tmp[cellid %in% sample.anno[qc.rna==T &qc.dna==T &!sex=="GG" ,cellid] ]
# acc
tmp<-tmp[cellid %in% sample.anno[qc.rna==T &qc.dna==T & qc.acc.extra==T,cellid] ]

#pdf(paste0("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/random.set/20190812.met.level.icm.te.",ndrtypes,"_",genetypes,".w.random_stages.pdf"),width=4,height=2.5)
pdf(paste0("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/random.set/20190921.acc.level.icm.te.w.random_stages.pdf"),width=10,height=4)

pdf(paste0("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/random.set/20190921.met.level.icm.te.w.random_stages.pdf"),width=10,height=4)
p1<-ggplot(tmp,aes(x=window_center, y=epi.level*100, fill=stage, color=stage ,group=stage)) +
  stat_summary(data=tmp[anno=="CRE"],aes(fill=stage, color=stage), fun.data=mean_sd, geom="smooth", alpha=0.2, size=0.5) + #aes(fill=stage, color=stage),
  stat_summary(data=tmp[anno=="random"],aes(fill=stage, color=stage), fun.data=mean_sd, geom="smooth", alpha=0.2, size=0.2,linetype="dashed") +
  scale_color_manual(values = c(  zygote="#7FC97F",`2cell`= "#BEAED4", `4cell`= "#FDC086",L4cell= "#808000",
                                  `8cell`= "#386CB0",`16cell`= "#650136", ICM= "#E64B35CC",TE= "#4DBBD5CC"))+ 
  scale_x_continuous(breaks=c(-2925,0,2925),labels=c("-3 kb","NDR","+3 kb"))+
  geom_vline(xintercept=0, linetype="dashed", color="grey", size=0.2) +
  labs(y="Met level (%)",x="") +
  coord_cartesian(ylim =c(0,50))+ 
  #labs(y="Acc level (%)",x="") +
  facet_grid(ndr.types~sig.type+gene.type)+
  #coord_cartesian(ylim =c(0,80))+ 
  theme_classic()+
  theme(
    axis.text = element_text(size=12),
    axis.title = element_text(size=14),
    legend.justification = "center",
    legend.key.width = unit(2.0,"line"),
    legend.title = element_blank(),
    legend.text = element_text(size=12),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(size=12),
    strip.background = element_blank()
  )
print(p1)
dev.off()
#  }  
#}


library(data.table)
library(ggplot2)
library(dplyr)


acc<-fread("~/Desktop/NGS_postdoc/data/DNA.data/201903.parsed.epi/20190529.sc.acc.NDRs.hiAccRNA.w.N.Nacc.tsv.gz")
met<-fread("~/Desktop/NGS_postdoc/data/DNA.data/201903.parsed.epi/20190529.sc.met.NDRs.hiAccRNA.w.N.Nacc.tsv.gz")

accrna.cor<-rbind(fread("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/20190810.distal.NDR_accrnaTPM.wtd.cor.final.tsv") %>% .[,anno:="Distal"] %>% .[!sig.type=="not.sig",.(id,anno,symbol,gene.type,sig.type)]#sig.type=="pos.sig"
                  ,fread("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/20190810.proximal.NDR_accrnaTPM.wtd.cor.final.tsv") %>% .[,anno:="Promoter"] %>% .[!sig.type=="not.sig",.(id,anno,symbol,gene.type,sig.type)] #sig.type=="pos.sig"
)

#accrna.cor<-rbind(fread("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/20190810.distal.NDR_accrnaTPM.wtd.cor.final.tsv") %>% .[,anno:="Distal"] %>% .[sig.type=="pos.sig",.(id,anno,symbol,gene.type,sig.type)]
#                  ,fread("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/20190810.proximal.NDR_accrnaTPM.wtd.cor.final.tsv") %>% .[,anno:="Promoter"] %>% .[sig.type=="pos.sig",.(id,anno,symbol,gene.type,sig.type)]
#)

sample.anno<-fread("~/Desktop/NGS_postdoc/data/DNA.data/20190524.updated.cell.metadata.tsv")
sample.anno[,gch_mean_rate:=gch_mean_rate %>% sub("%","",.) %>% as.numeric()]
sample.anno[,wcg_mean_rate:=wcg_mean_rate %>% sub("%","",.) %>% as.numeric()]

acc<-merge(acc,sample.anno[qc.rna==T &qc.dna==T &qc.acc.extra==T,.(sample=DNA.info,cellid,stage)],by="sample")
met<-merge(met,sample.anno[qc.rna==T &qc.dna==T & !sex=="GG",.(sample=DNA.info,cellid,stage)],by="sample")

accrna.cor.acc<-merge(accrna.cor,acc,by="id")
accrna.cor.met<-merge(accrna.cor,met,by="id")
tmp

accrna.cor

to.plot.mean<-accrna.cor.acc[weight>=3,.(mean(mean_rate)/100),.(cellid,stage,anno,gene.type,sig.type)] %>% merge(sample.anno[,.(cellid,bg=wcg_mean_rate/100)],by="cellid")

to.plot.mean$stage<-factor(to.plot.mean$stage,levels = c("zygote","2cell","4cell","L4cell","8cell","16cell","ICM","TE"))
to.plot.mean$anno<-factor(to.plot.mean$anno,levels = c("Promoter","Distal"))
to.plot.mean$gene.type<-factor(to.plot.mean$gene.type,levels = c("ICM.high","TE.high"))
to.plot.mean$sig.type<-factor(to.plot.mean$sig.type,levels = c("pos.sig","neg.sig"))
tmp<-dcast(to.plot.mean,cellid+stage+anno+sig.type~gene.type,value.var = "V1")

pdf("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/20190821.CREs.acc.te.minus.icm.pdf",width=10,height=4)
ggplot(tmp,aes(stage,(TE.high-ICM.high)*100,color=stage))+
  geom_jitter(alpha=0.8,size=0.5)+
  geom_violin(fill=NA,size=1,alpha=0.8)+
  geom_boxplot(width=0.1,fill=NA,color="black",alpha=0.5,outlier.colour = NA)+
  facet_grid(anno~sig.type)+
  geom_hline(yintercept = 0,linetype="dashed",color="black")+
  scale_color_manual(values = c(  zygote="#7FC97F",`2cell`= "#BEAED4", `4cell`= "#FDC086",L4cell= "#808000",
                                  `8cell`= "#386CB0",`16cell`= "#650136", ICM= "#E64B35CC",TE= "#4DBBD5CC"))+ # F0027F
  #geom_abline(slope = 1,linetype="dashed",color="grey")+
  theme_bw()+
  labs(x="",y="Differential Acc level (%)\nTE.CREs-ICM.CREs")+
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size=12),
    axis.text.x = element_text(size=12,angle=45,hjust=1),
    axis.title = element_text(size=14),
    legend.justification = "center",
    legend.title = element_blank(),
    legend.text = element_text(size=12),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank()
  ) 
dev.off()

to.plot.mean<-accrna.cor.met[weight>=3,.(mean(mean_rate)/100),.(cellid,stage,anno,gene.type,sig.type)] %>% merge(sample.anno[,.(cellid,bg=wcg_mean_rate/100)],by="cellid")

to.plot.mean$stage<-factor(to.plot.mean$stage,levels = c("zygote","2cell","4cell","L4cell","8cell","16cell","ICM","TE"))
to.plot.mean$anno<-factor(to.plot.mean$anno,levels = c("Promoter","Distal"))
to.plot.mean$gene.type<-factor(to.plot.mean$gene.type,levels = c("ICM.high","TE.high"))
to.plot.mean$sig.type<-factor(to.plot.mean$sig.type,levels = c("pos.sig","neg.sig"))
tmp<-dcast(to.plot.mean,cellid+stage+anno+sig.type~gene.type,value.var = "V1")

pdf("~/Desktop/NGS_postdoc/data/te.icm/icm.te.degs.cor/beforediff/20190821.CREs.met.te.minus.icm.pdf",width=10,height=4)
ggplot(tmp,aes(stage,(TE.high-ICM.high)*100,color=stage))+
  geom_jitter(alpha=0.8,size=0.5)+
  geom_violin(fill=NA,size=1,alpha=0.8)+
  geom_boxplot(width=0.1,fill=NA,color="black",alpha=0.5,outlier.colour = NA)+
  facet_grid(anno~sig.type)+
  geom_hline(yintercept = 0,linetype="dashed",color="black")+
  scale_color_manual(values = c(  zygote="#7FC97F",`2cell`= "#BEAED4", `4cell`= "#FDC086",L4cell= "#808000",
                                  `8cell`= "#386CB0",`16cell`= "#650136", ICM= "#E64B35CC",TE= "#4DBBD5CC"))+ # F0027F
  #geom_abline(slope = 1,linetype="dashed",color="grey")+
  theme_bw()+
  labs(x="",y="Differential Met level (%)\nTE.CREs-ICM.CREs")+
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size=12),
    axis.text.x = element_text(size=12,angle=45,hjust=1),
    axis.title = element_text(size=14),
    legend.justification = "center",
    legend.title = element_blank(),
    legend.text = element_text(size=12),
    # panel.border = element_blank(),
    # panel.grid.major = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank()
  ) 
dev.off()


