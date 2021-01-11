# figure 3c-h

##### plot allelic epigenomes vs total RNA
# 20190715

# run in server 
# /home/wangyang/project/duozuxue/DNA.part/local_correlation/allelic.cor


# plot in local

##########################################################
# plot whole genic local correlations 
##########################################################

# merge tss genebody and tts
opts<-list()
opts$up <- 2500
opts$window <- 500
opts$slide <- 100
opts$down <- 2500
tmp <- seq(from=0-opts$up, to=0+opts$down-opts$window, by=opts$slide)
tss.order <- data.table(window_center=tmp+(opts$window/2), rel_start=tmp, rel_end=tmp+opts$window)
tss.order<-tss.order[window_center<=50]
tss.order[,reord:=.I]
tss.order[,anno:="tss2500"]

opts$up <- 0
opts$window <- 10
opts$slide <- 2
opts$down <- 100
tmp <- seq(from=0-opts$up, to=0+opts$down-opts$window, by=opts$slide)
genebody.order <- data.table(window_center=tmp+(opts$window/2), rel_start=tmp, rel_end=tmp+opts$window)
genebody.order[,reord:=.I+24]
genebody.order[,anno:="genebody"]


opts$up <- 2500
opts$window <- 500
opts$slide <- 100
opts$down <- 2500
tmp <- seq(from=0-opts$up, to=0+opts$down-opts$window, by=opts$slide)
tts.order <- data.table(window_center=tmp+(opts$window/2), rel_start=tmp, rel_end=tmp+opts$window)
tts.order<-tts.order[window_center>= -75]
tts.order[,reord:=.I+70]
tts.order[,anno:="tts2500"]

new.order<-rbind(tss.order,genebody.order,tts.order) %>% .[,.(window_center,reord,anno)]

for (epi.type in c("acc", "met")){
  epi.type="met"  
files<-list.files("~/Desktop/NGS_postdoc/data/correlations/local.correlation/allelic/20190715.allelic.epi.vs.all.RNA",pattern = paste0("20190715.cor.local.",epi.type,"rna."),full.names = T)
cor<-lapply(files,function(f)fread(f) %>% .[,stage:=basename(f) %>% sub( paste0("20190715.cor.local.",epi.type,"rna."),"",.) %>% sub(".tsv","",.)]) %>% rbindlist()

cor$allelic.type<-do.call(rbind,strsplit(cor$stage,"[.]"))[,2]
cor$anno<-do.call(rbind,strsplit(cor$stage,"[.]"))[,3]
cor$stage<-do.call(rbind,strsplit(cor$stage,"[.]"))[,1]

unique(cor$anno)
unique(cor$stage)

cor$stage<-factor(cor$stage,levels = c("zygote","2cell","4cell","L4cell","8cell","16cell","ICM","TE"))
cor$anno<-factor(cor$anno,levels=c("tss2500","genebody","tts2500"))
# merge tss genebody and tts
cor[anno=="tss2500"]


cor<-merge(cor,new.order,by=c("window_center","anno"))

sample.anno<-fread("~/Desktop/NGS_postdoc/data/DNA.data/20190524.updated.cell.metadata.tsv")

# select cells to plot
if (epi.type=="met"){
  cor<-merge(cor,sample.anno[qc.rna==T & qc.dna==T & !sex=="GG",.(cellid,qc.acc.extra)],by="cellid")
} else {
  cor<-merge(cor,sample.anno[qc.rna==T & qc.dna==T &qc.acc.extra==T,.(cellid,qc.acc.extra)],by="cellid")
}

xlen<-cor[,length(unique(reord))]
ylen<-cor[,max(r)-min(r)]

pdf(paste0("~/Desktop/NGS_postdoc/data/correlations/local.correlation/plot/20190916.cor.local.allelic.",epi.type,".vs.all.rna.gene.pdf"),width=10,height=4) 
p1<-ggplot(cor,aes(x=reord, y=r, fill=stage, color=stage)) +
  stat_summary(aes(fill=stage, color=stage), fun.data=mean_sd, geom="smooth", alpha=0.5, size=0.75) + #aes(fill=stage, color=stage),
  scale_color_manual(values = c(  zygote="#7FC97F",`2cell`= "#BEAED4", `4cell`= "#FDC086",L4cell= "#808000",
                                  `8cell`= "#386CB0",`16cell`= "#650136", ICM= "#E64B35CC",TE= "#4DBBD5CC"))+
  geom_hline(yintercept=0, linetype="dashed", color="black", size=0.5) +
  geom_vline(xintercept=24, color="grey40", size=0.5, alpha=0.5,linetype="dashed") +
  geom_vline(xintercept=70,  color="grey40", size=0.5,alpha=0.5,linetype="dashed") +
  facet_grid(.~allelic.type)+
  ylab(paste0(epi.type,"/RNA Correlation")) + xlab("") +
  scale_x_continuous(breaks = c(1,24,47,70,94),labels = c("-2.5 kb","TSS","Genebody","TES","2.5kb"))+
  theme_bw()+
  theme(
    # axis.text.x = element_text(angle = 45,hjust=1),
    #   axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=rel(1.1), colour="black"),
    axis.title.y = element_text(size=rel(1.2)),
    legend.key = element_blank(),
    legend.position = "right",
    legend.justification = "center",
    legend.text = element_text(size=rel(1.2)),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size=rel(1.2))
  )+coord_fixed(ratio = xlen/ylen)
print (p1)
dev.off()
}


###############################################
# plot levels of allelic acc met
###############################################
for (epi.type in c("acc", "met")){
 # epi.type="acc"
  files<-list.files("~/Desktop/NGS_postdoc/data/correlations/local.correlation/allelic/20190715.allelic.epi.vs.all.RNA",pattern = paste0("20190715.local.",epi.type,".level."),full.names = T)
  cor<-lapply(files,function(f)fread(f) %>% .[,stage:=basename(f) %>% sub( paste0("20190715.local.",epi.type,".level."),"",.) %>% sub(".tsv","",.)]) %>% rbindlist()
  
  cor$allelic.type<-do.call(rbind,strsplit(cor$stage,"[.]"))[,2]
  cor$anno<-do.call(rbind,strsplit(cor$stage,"[.]"))[,3]
  cor$stage<-do.call(rbind,strsplit(cor$stage,"[.]"))[,1]
  
  unique(cor$anno)
  unique(cor$stage)
  
  cor$stage<-factor(cor$stage,levels = c("zygote","2cell","4cell","L4cell","8cell","16cell","ICM","TE"))
  cor$anno<-factor(cor$anno,levels=c("tss2500","genebody","tts2500"))
  # merge tss genebody and tts
  cor[anno=="tss2500"]
  
  
  cor<-merge(cor,new.order,by=c("window_center","anno"))
  
  sample.anno<-fread("~/Desktop/NGS_postdoc/data/DNA.data/20190524.updated.cell.metadata.tsv")
  
  # select cells to plot
  if (epi.type=="met"){
    cor<-merge(cor,sample.anno[qc.rna==T & qc.dna==T & !sex=="GG",.(cellid,qc.acc.extra)],by="cellid")
  } else {
    cor<-merge(cor,sample.anno[qc.rna==T & qc.dna==T &qc.acc.extra==T,.(cellid,qc.acc.extra)],by="cellid")
  }
  xlen<-cor[,length(unique(reord))]
  ylen<-cor[,max(epi.level*100)-min(epi.level*100)]
  pdf(paste0("~/Desktop/NGS_postdoc/data/correlations/local.correlation/plot/20190916.local.allelic.",epi.type,".level.pdf"),width=10,height=4) 
  p1<-ggplot(cor,aes(x=reord, y=epi.level*100, fill=stage, color=stage)) +
    stat_summary(aes(fill=stage, color=stage), fun.data=mean_sd, geom="smooth", alpha=0.5, size=0.75) + #aes(fill=stage, color=stage),
    scale_color_manual(values = c(  zygote="#7FC97F",`2cell`= "#BEAED4", `4cell`= "#FDC086",L4cell= "#808000",
                                    `8cell`= "#386CB0",`16cell`= "#650136", ICM= "#E64B35CC",TE= "#4DBBD5CC"))+
    geom_vline(xintercept=24, color="grey40", size=0.5, alpha=0.5,linetype="dashed") +
    geom_vline(xintercept=70,  color="grey40", size=0.5,alpha=0.5,linetype="dashed") +
    facet_grid(.~allelic.type)+
    ylab(paste0(epi.type," level %")) + xlab("") +
    scale_x_continuous(breaks = c(1,24,47,70,94),labels = c("-2.5 kb","TSS","Genebody","TES","2.5kb"))+
    theme_bw()+
    theme(
      #axis.text.x = element_text(angle = 45,hjust=1),
      #   axis.ticks.x = element_blank(),
      axis.text.y = element_text(size=rel(1.1), colour="black"),
      axis.title.y = element_text(size=rel(1.2)),
      legend.key = element_blank(),
      legend.position = "right",
      legend.justification = "center",
      legend.text = element_text(size=rel(1.2)),
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size=rel(1.2))
    )+coord_fixed(ratio = xlen/ylen)
  print (p1)
  dev.off()
}


##############  
## allelic epi vs total rna rm maternal  20950 


# local allelic correlations 
# run in server 
# /home/wangyang/project/duozuxue/DNA.part/local_correlation/allelic.cor/20190915.allelic.epi.vs.all.rna.rm.maternal.gene
#20190916.allelic.epi.vs.all.rna.only.maternal.gene


for (epi.type in c("acc", "met")){
  epi.type="met"
  files<-list.files("~/Desktop/NGS_postdoc/data/correlations/local.correlation/allelic/allelic.cor/20190915.allelic.epi.vs.all.rna.rm.maternal.gene",pattern = paste0("20190915.cor.local.",epi.type,"rna."),full.names = T)
  cor<-lapply(files,function(f)fread(f) %>% .[,stage:=basename(f) %>% sub( paste0("20190915.cor.local.",epi.type,"rna."),"",.) %>% sub(".tsv","",.)]) %>% rbindlist()
  
cor$allelic.type<-do.call(rbind,strsplit(cor$stage,"[.]"))[,2]
cor$anno<-do.call(rbind,strsplit(cor$stage,"[.]"))[,3]
cor$stage<-do.call(rbind,strsplit(cor$stage,"[.]"))[,1]

unique(cor$anno)
#cor[,anno:=anno%>% sub("zygote|2cell|4cell|8cell|L4cell|8cell|16cell|ICM|TE","",.)]
#cor[,stage:=stage%>% sub("genebody|tts2500","",.)]
unique(cor$stage)

cor$stage<-factor(cor$stage,levels = c("zygote","2cell","4cell","L4cell","8cell","16cell","ICM","TE"))
cor$anno<-factor(cor$anno,levels=c("tss2500","genebody","tts2500"))
# merge tss genebody and tts
cor[anno=="tss2500"]

cor<-merge(cor,new.order,by=c("window_center","anno"))
sample.anno<-fread("~/Desktop/NGS_postdoc/data/DNA.data/20190524.updated.cell.metadata.tsv")

# select cells to plot
if (epi.type=="met"){
  cor<-merge(cor,sample.anno[qc.rna==T & qc.dna==T & !sex=="GG",.(cellid,qc.acc.extra)],by="cellid")
} else {
  cor<-merge(cor,sample.anno[qc.rna==T & qc.dna==T &qc.acc.extra==T,.(cellid,qc.acc.extra)],by="cellid")
}

xlen<-cor[,length(unique(reord))]
ylen<-cor[,max(r)-min(r)]

library(ggpubr) # mean_sd
xlen<-cor[,length(unique(reord))]
ylen<-cor[,max(r)-min(r)]
ylen<- 0.6
pdf(paste0("~/Desktop/NGS_postdoc/data/correlations/local.correlation/plot/20190916.cor.local.allelic.",epi.type,".vs.all.rna.gene.rm.maternal.gene.pdf"),width=4,height=10) 
p1<-ggplot(cor,aes(x=reord, y=r, fill=stage, color=stage)) +
  stat_summary(aes(fill=cellid, color=stage), fun.data=mean_sd, geom="smooth", alpha=0.5, size=0.75) + #aes(fill=stage, color=stage),
  scale_color_manual(values = c(  zygote="#7FC97F",`2cell`= "#BEAED4", `4cell`= "#FDC086",L4cell= "#808000",
                                  `8cell`= "#386CB0",`16cell`= "#650136", ICM= "#E64B35CC",TE= "#4DBBD5CC"))+
  stat_summary(aes(fill=stage, color=stage), fun.data=mean_sd, geom="smooth", alpha=0.5, size=0.75,color="black") + #aes(fill=stage, color=stage),
  geom_hline(yintercept=0, linetype="dashed", color="black", size=0.5) +
  geom_vline(xintercept=24, color="grey40", size=0.5, alpha=0.5,linetype="dashed") +
  geom_vline(xintercept=70,  color="grey40", size=0.5,alpha=0.5,linetype="dashed") +
  facet_grid(stage~allelic.type)+
  guides(fill=F)+
  #ylab("Met/RNA Correlation") + xlab("") +
  ylab(paste0(epi.type,"/RNA Correlation")) + xlab("") +
  scale_x_continuous(breaks = c(1,24,47,70,94),labels = c("-2.5 kb","TSS","Genebody","TES","+2.5 kb"))+
  theme_bw()+
  theme(
     axis.text.x = element_text(angle = 45,hjust=1),
    #   axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=rel(1.1), colour="black"),
    axis.title.y = element_text(size=rel(1.2)),
    legend.key = element_blank(),
    legend.position = "right",
    legend.justification = "center",
    legend.text = element_text(size=rel(1.2)),
    legend.title = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size=rel(1.2))
  )+coord_fixed(ratio = xlen/ylen,ylim = c(-0.3,0.3))
print (p1)
dev.off()
}

###############################################
# plot levels of acc met
############################################### 20190915.local.met.level.
for (epi.type in c("acc", "met")){
  files<-list.files("~/Desktop/NGS_postdoc/data/correlations/local.correlation/allelic/allelic.cor/20190915.allelic.epi.vs.all.rna.rm.maternal.gene",pattern = paste0("20190915.local.",epi.type,".level."),full.names = T)
  cor<-lapply(files,function(f)fread(f) %>% .[,stage:=basename(f) %>% sub( paste0("20190915.local.",epi.type,".level."),"",.) %>% sub(".tsv","",.)]) %>% rbindlist()
  
  cor$allelic.type<-do.call(rbind,strsplit(cor$stage,"[.]"))[,2]
  cor$anno<-do.call(rbind,strsplit(cor$stage,"[.]"))[,3]
  cor$stage<-do.call(rbind,strsplit(cor$stage,"[.]"))[,1]
  
  unique(cor$anno)
  #cor[,anno:=anno%>% sub("zygote|2cell|4cell|8cell|L4cell|8cell|16cell|ICM|TE","",.)]
  #cor[,stage:=stage%>% sub("genebody|tts2500","",.)]
  unique(cor$stage)
  
  cor$stage<-factor(cor$stage,levels = c("zygote","2cell","4cell","L4cell","8cell","16cell","ICM","TE"))
  cor$anno<-factor(cor$anno,levels=c("tss2500","genebody","tts2500"))
  # merge tss genebody and tts
  cor[anno=="tss2500"]
  
  cor<-merge(cor,new.order,by=c("window_center","anno"))
  sample.anno<-fread("~/Desktop/NGS_postdoc/data/DNA.data/20190524.updated.cell.metadata.tsv")
  
  # select cells to plot
  if (epi.type=="met"){
    cor<-merge(cor,sample.anno[qc.rna==T & qc.dna==T & !sex=="GG",.(cellid,qc.acc.extra)],by="cellid")
  } else {
    cor<-merge(cor,sample.anno[qc.rna==T & qc.dna==T &qc.acc.extra==T,.(cellid,qc.acc.extra)],by="cellid")
  }
  
  library(ggpubr) # mean_sd
  xlen<-cor[,length(unique(reord))]
  ylen<-cor[,max(epi.level*100)-min(epi.level*100)]
  ylen= 80
  pdf(paste0("~/Desktop/NGS_postdoc/data/correlations/local.correlation/plot/20190916.local.allelic.",epi.type,".level.rm.maternal.gene.pdf"),width=4,height=10) 
  p1<-ggplot(cor,aes(x=reord, y=epi.level*100, fill=stage, color=stage)) +
    stat_summary(aes(fill=cellid, color=stage), fun.data=mean_sd, geom="smooth", alpha=0.5, size=0.75) + #aes(fill=stage, color=stage),
    scale_color_manual(values = c(  zygote="#7FC97F",`2cell`= "#BEAED4", `4cell`= "#FDC086",L4cell= "#808000",
                                    `8cell`= "#386CB0",`16cell`= "#650136", ICM= "#E64B35CC",TE= "#4DBBD5CC"))+
    stat_summary(aes(fill=stage, color=stage), fun.data=mean_sd, geom="smooth", alpha=0.5, size=0.75,color="black") + #aes(fill=stage, color=stage),
 #   geom_hline(yintercept=0, linetype="dashed", color="black", size=0.5) +
    geom_vline(xintercept=24, color="grey40", size=0.5, alpha=0.5,linetype="dashed") +
    geom_vline(xintercept=70,  color="grey40", size=0.5,alpha=0.5,linetype="dashed") +
    facet_grid(stage~allelic.type)+
    guides(fill=F)+
    ylab(paste0(epi.type," level %")) + xlab("") +
    scale_x_continuous(breaks = c(1,24,47,70,94),labels = c("-2.5 kb","TSS","Genebody","TES","+2.5 kb"))+
    theme_bw()+
    theme(
      axis.text.x = element_text(angle = 45,hjust=1),
      #   axis.ticks.x = element_blank(),
      axis.text.y = element_text(size=rel(1.1), colour="black"),
      axis.title.y = element_text(size=rel(1.2)),
      legend.key = element_blank(),
      legend.position = "right",
      legend.justification = "center",
      legend.text = element_text(size=rel(1.2)),
      legend.title = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size=rel(1.2))
    )+coord_fixed(ratio = xlen/ylen,ylim=c(0,80))
  print (p1)
  dev.off()
}



##############  
## allelic epi vs total rna only maternal (average expr>=1 at zygote )   12939


# local allelic correlations 
# run in server 
# /home/wangyang/project/duozuxue/DNA.part/local_correlation/allelic.cor/20190916.allelic.epi.vs.all.rna.only.maternal.gene


for (epi.type in c("acc", "met")){
  files<-list.files("~/Desktop/NGS_postdoc/data/correlations/local.correlation/allelic/allelic.cor/20190916.allelic.epi.vs.all.rna.only.maternal.gene",pattern = paste0("20190916.cor.local.",epi.type,"rna."),full.names = T)
  cor<-lapply(files,function(f)fread(f) %>% .[,stage:=basename(f) %>% sub( paste0("20190916.cor.local.",epi.type,"rna."),"",.) %>% sub(".tsv","",.)]) %>% rbindlist()
  
  cor$allelic.type<-do.call(rbind,strsplit(cor$stage,"[.]"))[,2]
  cor$anno<-do.call(rbind,strsplit(cor$stage,"[.]"))[,3]
  cor$stage<-do.call(rbind,strsplit(cor$stage,"[.]"))[,1]
  
  unique(cor$anno)
  #cor[,anno:=anno%>% sub("zygote|2cell|4cell|8cell|L4cell|8cell|16cell|ICM|TE","",.)]
  #cor[,stage:=stage%>% sub("genebody|tts2500","",.)]
  unique(cor$stage)
  
  cor$stage<-factor(cor$stage,levels = c("zygote","2cell","4cell","L4cell","8cell","16cell","ICM","TE"))
  cor$anno<-factor(cor$anno,levels=c("tss2500","genebody","tts2500"))
  # merge tss genebody and tts
  cor[anno=="tss2500"]
  
  cor<-merge(cor,new.order,by=c("window_center","anno"))
  sample.anno<-fread("~/Desktop/NGS_postdoc/data/DNA.data/20190524.updated.cell.metadata.tsv")
  
  # select cells to plot
  if (epi.type=="met"){
    cor<-merge(cor,sample.anno[qc.rna==T & qc.dna==T & !sex=="GG",.(cellid,qc.acc.extra)],by="cellid")
  } else {
    cor<-merge(cor,sample.anno[qc.rna==T & qc.dna==T &qc.acc.extra==T,.(cellid,qc.acc.extra)],by="cellid")
  }
  
  xlen<-cor[,length(unique(reord))]
  ylen<-cor[,max(r)-min(r)]
  
  library(ggpubr) # mean_sd
  
  pdf(paste0("~/Desktop/NGS_postdoc/data/correlations/local.correlation/plot/20190916.cor.local.allelic.",epi.type,".vs.all.rna.gene.only.maternal.gene.pdf"),width=4,height=10) 
  p1<-ggplot(cor,aes(x=reord, y=r, fill=stage, color=stage)) +
    stat_summary(aes(fill=cellid, color=stage), fun.data=mean_sd, geom="smooth", alpha=0.5, size=0.75) + #aes(fill=stage, color=stage),
    scale_color_manual(values = c(  zygote="#7FC97F",`2cell`= "#BEAED4", `4cell`= "#FDC086",L4cell= "#808000",
                                    `8cell`= "#386CB0",`16cell`= "#650136", ICM= "#E64B35CC",TE= "#4DBBD5CC"))+
    stat_summary(aes(fill=stage, color=stage), fun.data=mean_sd, geom="smooth", alpha=0.5, size=0.75,color="black") + #aes(fill=stage, color=stage),
    geom_hline(yintercept=0, linetype="dashed", color="black", size=0.5) +
    geom_vline(xintercept=24, color="grey40", size=0.5, alpha=0.5,linetype="dashed") +
    geom_vline(xintercept=70,  color="grey40", size=0.5,alpha=0.5,linetype="dashed") +
    facet_grid(stage~allelic.type)+
    guides(fill=F)+
    #ylab("Met/RNA Correlation") + xlab("") +
    ylab(paste0(epi.type,"/RNA Correlation")) + xlab("") +
    scale_x_continuous(breaks = c(1,24,47,70,94),labels = c("-2.5 kb","TSS","Genebody","TES","+2.5 kb"))+
    theme_bw()+
    theme(
      axis.text.x = element_text(angle = 45,hjust=1),
      #   axis.ticks.x = element_blank(),
      axis.text.y = element_text(size=rel(1.1), colour="black"),
      axis.title.y = element_text(size=rel(1.2)),
      legend.key = element_blank(),
      legend.position = "right",
      legend.justification = "center",
      legend.text = element_text(size=rel(1.2)),
      legend.title = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size=rel(1.2))
    )+coord_fixed(ratio = xlen/ylen, ylim = c(-0.3,0.3))
  print (p1)
  dev.off()
}

###############################################
# plot levels of acc met
############################################### 20190915.local.met.level.
for (epi.type in c("acc", "met")){
  files<-list.files("~/Desktop/NGS_postdoc/data/correlations/local.correlation/allelic/allelic.cor/20190916.allelic.epi.vs.all.rna.only.maternal.gene",pattern = paste0("20190916.local.",epi.type,".level."),full.names = T)
  cor<-lapply(files,function(f)fread(f) %>% .[,stage:=basename(f) %>% sub( paste0("20190916.local.",epi.type,".level."),"",.) %>% sub(".tsv","",.)]) %>% rbindlist()
  
  cor$allelic.type<-do.call(rbind,strsplit(cor$stage,"[.]"))[,2]
  cor$anno<-do.call(rbind,strsplit(cor$stage,"[.]"))[,3]
  cor$stage<-do.call(rbind,strsplit(cor$stage,"[.]"))[,1]
  
  unique(cor$anno)
  #cor[,anno:=anno%>% sub("zygote|2cell|4cell|8cell|L4cell|8cell|16cell|ICM|TE","",.)]
  #cor[,stage:=stage%>% sub("genebody|tts2500","",.)]
  unique(cor$stage)
  
  cor$stage<-factor(cor$stage,levels = c("zygote","2cell","4cell","L4cell","8cell","16cell","ICM","TE"))
  cor$anno<-factor(cor$anno,levels=c("tss2500","genebody","tts2500"))
  # merge tss genebody and tts
  cor[anno=="tss2500"]
  
  cor<-merge(cor,new.order,by=c("window_center","anno"))
  sample.anno<-fread("~/Desktop/NGS_postdoc/data/DNA.data/20190524.updated.cell.metadata.tsv")
  
  # select cells to plot
  if (epi.type=="met"){
    cor<-merge(cor,sample.anno[qc.rna==T & qc.dna==T & !sex=="GG",.(cellid,qc.acc.extra)],by="cellid")
  } else {
    cor<-merge(cor,sample.anno[qc.rna==T & qc.dna==T &qc.acc.extra==T,.(cellid,qc.acc.extra)],by="cellid")
  }
  
  library(ggpubr) # mean_sd
  xlen<-cor[,length(unique(reord))]
  ylen<-cor[,max(epi.level*100)-min(epi.level*100)]
  
  pdf(paste0("~/Desktop/NGS_postdoc/data/correlations/local.correlation/plot/20190916.local.allelic.",epi.type,".level.only.maternal.gene.pdf"),width=4,height=10) 
  p1<-ggplot(cor,aes(x=reord, y=epi.level*100, fill=stage, color=stage)) +
    stat_summary(aes(fill=cellid, color=stage), fun.data=mean_sd, geom="smooth", alpha=0.5, size=0.75) + #aes(fill=stage, color=stage),
    scale_color_manual(values = c(  zygote="#7FC97F",`2cell`= "#BEAED4", `4cell`= "#FDC086",L4cell= "#808000",
                                    `8cell`= "#386CB0",`16cell`= "#650136", ICM= "#E64B35CC",TE= "#4DBBD5CC"))+
    stat_summary(aes(fill=stage, color=stage), fun.data=mean_sd, geom="smooth", alpha=0.5, size=0.75,color="black") + #aes(fill=stage, color=stage),
    #   geom_hline(yintercept=0, linetype="dashed", color="black", size=0.5) +
    geom_vline(xintercept=24, color="grey40", size=0.5, alpha=0.5,linetype="dashed") +
    geom_vline(xintercept=70,  color="grey40", size=0.5,alpha=0.5,linetype="dashed") +
    facet_grid(stage~allelic.type)+
    guides(fill=F)+
    ylab(paste0(epi.type," level %")) + xlab("") +
    scale_x_continuous(breaks = c(1,24,47,70,94),labels = c("-2.5 kb","TSS","Genebody","TES","+2.5 kb"))+
    theme_bw()+
    theme(
      axis.text.x = element_text(angle = 45,hjust=1),
      #   axis.ticks.x = element_blank(),
      axis.text.y = element_text(size=rel(1.1), colour="black"),
      axis.title.y = element_text(size=rel(1.2)),
      legend.key = element_blank(),
      legend.position = "right",
      legend.justification = "center",
      legend.text = element_text(size=rel(1.2)),
      legend.title = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size=rel(1.2))
    )+coord_fixed(ratio = xlen/ylen,ylim = c(0,80))
  print (p1)
  dev.off()
}


#################################################20200519

# statistical analysis maternal and nonmaternal gene genebody correlations
epi.type="met"

files<-list.files("~/Desktop/NGS_postdoc/data/correlations/local.correlation/allelic/allelic.cor/20190915.allelic.epi.vs.all.rna.rm.maternal.gene",pattern = paste0("20190915.cor.local.",epi.type,"rna."),full.names = T)
cor.nonmaternal<-lapply(files,function(f)fread(f) %>% .[,stage:=basename(f) %>% sub( paste0("20190915.cor.local.",epi.type,"rna."),"",.) %>% sub(".tsv","",.)]) %>% rbindlist() %>% .[,gene.type:="nonmaternal"]

files<-list.files("~/Desktop/NGS_postdoc/data/correlations/local.correlation/allelic/allelic.cor/20190916.allelic.epi.vs.all.rna.only.maternal.gene",pattern = paste0("20190916.cor.local.",epi.type,"rna."),full.names = T)
cor.maternal<-lapply(files,function(f)fread(f) %>% .[,stage:=basename(f) %>% sub( paste0("20190916.cor.local.",epi.type,"rna."),"",.) %>% sub(".tsv","",.)]) %>% rbindlist()%>% .[,gene.type:="maternal"]

cor<-rbind(cor.nonmaternal,cor.maternal)


cor$allelic.type<-do.call(rbind,strsplit(cor$stage,"[.]"))[,2]
cor$anno<-do.call(rbind,strsplit(cor$stage,"[.]"))[,3]
cor$stage<-do.call(rbind,strsplit(cor$stage,"[.]"))[,1]

unique(cor$anno)
#cor[,anno:=anno%>% sub("zygote|2cell|4cell|8cell|L4cell|8cell|16cell|ICM|TE","",.)]
#cor[,stage:=stage%>% sub("genebody|tts2500","",.)]
unique(cor$stage)

cor$stage<-factor(cor$stage,levels = c("zygote","2cell","4cell","L4cell","8cell","16cell","ICM","TE"))
cor$anno<-factor(cor$anno,levels=c("tss2500","genebody","tts2500"))

sample.anno<-fread("~/Desktop/NGS_postdoc/data/DNA.data/20190524.updated.cell.metadata.tsv")

# select cells to plot
cor<-merge(cor,sample.anno[qc.rna==T & qc.dna==T & !sex=="GG",.(cellid,qc.acc.extra)],by="cellid")

length(unique(cor$cellid))

to.plot<-cor[,mean(r),.(cellid,stage,gene.type,allelic.type,anno)]

pdf(paste0("~/Desktop/NGS_postdoc/data/correlations/local.correlation/plot/20200519.local.allelic.genebody.met.vs.expr.p.value.box.pdf"),width=12,height=4) 
ggplot(to.plot[anno=="genebody"],aes(gene.type,V1,color=gene.type))+
  geom_boxplot(outlier.colour = NA)+
  labs(x="",y="Allelic Met vs Expr Correlation")+
  stat_compare_means(label="p.format",method = "t.test",size=3,ref.group = "maternal")+
  facet_grid(allelic.type~stage)+
  theme_bw()+
  theme(
    axis.text.x=element_blank(),
    strip.background = element_blank()
  )
dev.off()