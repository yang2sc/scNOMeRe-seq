#library(GenomicRanges) 
#library(HMMcopy)
library(data.table)
library(purrr)
library(furrr)
library(dplyr)
library(BSgenome.Mmusculus.UCSC.mm9)
library(weights)
library(plyr)
library(RColorBrewer)

################################################################################################################
# DNA data
################################################################################################################
# remove failed rna, dna cells
cell.anno<-fread("~/Desktop/NGS_postdoc/data/DNA.data/20190524.updated.cell.metadata.tsv")
cell.anno<-cell.anno[qc.rna==T|qc.dna==T] %>% .[stage %in% c(  "zygote","2cell", "4cell","L4cell","8cell", "16cell","ICM","TE")]

io<-list()
io$copy<-"~/Desktop/NGS_postdoc/data/DNA.data/hmmcopy.cnv/data/"

test.copy <-  dir(io$copy,
                  recursive = F, full = TRUE, pattern = ".corrected.csv")%>%
  setNames(., basename(.) %>% sub(".hmmcopy.win1000000.corrected.csv", "", .) ) %>%
  map2(names(.), ~fread(.)[, anno := .y]) %>% rbindlist #%>% .[,c(2:5,9)] %>%setnames(c("motif", "id", "start", "stop","q"))


test.copy.plot<-test.copy[!is.na(copy)]
test.copy.plot<-test.copy.plot[,cor.map:=cor.map *2]

#test.copy.plot[substring(anno,1,5)=="later",anno:=paste0("L",substring(anno,7,20))]

### remove X Y chromsome
test.copy.plot<-test.copy.plot[space %in% paste0("chr",c(1:19))]

test.copy.plot<-merge(test.copy.plot[,.(chr=space,start,end,cor.map,anno)],cell.anno[,.(anno=cellid,stage,DNA.info,embryo_id,sex,cnv_info,qc.rna,qc.dna)],by="anno")
test.copy.plot$chr<-factor(test.copy.plot$chr,levels=c(paste0("chr",c(1:19))))
test.copy.plot$stage<-factor(test.copy.plot$stage,levels = c(  "zygote","2cell", "4cell","L4cell","8cell", "16cell","ICM","TE"))
test.copy.plot<- setkey(test.copy.plot,stage,cnv_info,sex,chr,start,end)

parse.metadata<-unique(test.copy.plot[,.(chr,start,end)])
parse.metadata[,id:=paste0("bin_",rownames(parse.metadata))]
#fwrite(parse.metadata,"~/Desktop/CNV&RNA&DNAmethylation/20190207.1M.bin.cnv.metadata.tsv",sep="\t")
#parse.metadata<-fread("~/Desktop/CNV&RNA&DNAmethylation/20190207.1M.bin.cnv.metadata.tsv")
setkey(parse.metadata,chr,start,end)


test.copy.plot<-merge(test.copy.plot,parse.metadata,by=c("chr","start","end"))
test.copy.plot.dcast<-dcast(test.copy.plot,anno~id,value.var = "cor.map")
# add cells rna 
test.copy.plot.dcast<-merge(test.copy.plot.dcast,cell.anno[,.(anno=cellid)],by="anno",all=T)
test.copy.plot.dcast<-test.copy.plot.dcast%>% as.data.frame()%>% tibble::column_to_rownames("anno")
test.copy.plot.dcast[1:5,1:5]
test.copy.plot.dcast[is.na(test.copy.plot.dcast)]
test.copy.plot.dcast[test.copy.plot.dcast>3]<- 3
test.copy.plot.dcast[test.copy.plot.dcast<1]<- 1

colanno<-parse.metadata[,.(id,chr)]%>% as.data.frame()%>% tibble::column_to_rownames("id") 
cell.anno[cellid=="X4cell6_4"]
rowanno<-cell.anno[cellid %in% rownames(test.copy.plot.dcast),.(cellid,stage,embryo_id,cnv_info,sex,GCH_sites)]%>% setkey(stage,embryo_id,cellid)%>% .[,embryo_id:=NULL]%>% as.data.frame()%>% tibble::column_to_rownames("cellid")
#chr10	3002075	+	0	1	0	GCT  #cnv_info,sex
col.gaps<-as.data.frame(table(colanno$chr))
col.gaps$group<-1
class(col.gaps)

col.gaps<- ddply(col.gaps, "group",
                 transform, sep=cumsum(Freq))

row.gaps<-as.data.frame(table(rowanno$stage))
row.gaps$group<-1
class(row.gaps)
#library(plyr)
row.gaps<- ddply(row.gaps, "group",
                 transform, sep=cumsum(Freq))

#cols <- colorRampPalette(brewer.pal(9, "Set1"))
#stagecolors <- cols(length(unique(rowanno$stage)))
#names(stagecolors) <- unique(rowanno$stage)

cols <- colorRampPalette(brewer.pal(8, "Set2"))
cnvcolors <- cols(length(unique(rowanno$cnv_info)))
names(cnvcolors) <- unique(rowanno$cnv_info)

cols <- colorRampPalette(brewer.pal(8, "Set3"))
sexcolors <- cols(length(unique(rowanno$sex)))
names(sexcolors) <- unique(rowanno$sex)

cols <- colorRampPalette(brewer.pal(9,"Blues"))
sitescolors <- cols(length(seq(0,max(rowanno$GCH_sites),by=5000000)))
names(sitescolors) <- unique(seq(0,max(rowanno$GCH_sites),by=5000000))


newCols <- colorRampPalette(grDevices::rainbow(length(unique(colanno$chr))))
cololors <- newCols(length(unique(colanno$chr)))
names(cololors) <- unique(colanno$chr)
annocolor <- list(chr = cololors,stage = c( zygote="#7FC97F",`2cell`= "#BEAED4", `4cell`= "#FDC086",L4cell= "#FFFF99",
                                            `8cell`= "#386CB0",`16cell`= "#F0027F", ICM= "#E64B35CC",TE= "#4DBBD5CC"),cnv_info=cnvcolors,sex=sexcolors,GCH_sites=sitescolors)



pdf(paste0("~/Desktop/NGS_postdoc/data/CNV/20191010.Genome_CNV.hmmcopy.pdf"),width = 15,height=15)
pheatmap::pheatmap(test.copy.plot.dcast[rownames(rowanno),rownames(colanno)],cluster_rows = F,cluster_cols = F,
                   cellwidth=0.25,cellheight = 4,
                   # main = i,
                   # col=rev(cols),
                   color = colorRampPalette(c("royalblue","white","orangered"))(100),
                   border_color = NA,
                   gaps_col=col.gaps$sep,gaps_row = row.gaps$sep,
                   annotation_col=colanno,
                   annotation_row = rowanno,
                   annotation_colors = annocolor,
                   # annotation_names_col=F,
                   show_rownames=F,show_colnames=F)

dev.off()

################################################################################################################
# RNA data
################################################################################################################
sample.anno<-fread("~/Desktop/NGS_postdoc/data/DNA.data/20190524.updated.cell.metadata.tsv") 
sample.anno[,unique(stage)]
sample.anno$stage<-factor(sample.anno$stage,levels = c(  "zygote","2cell", "4cell","L4cell","8cell", "16cell","ICM","TE"))

gene.anno<-fread("~/Desktop/NGS_postdoc/data/DNA.data/ref/wy.gencode.gene.metadata.tsv") %>% .[,.(chr,start,end,gene=ens_ID,symbol)] %>% setkey(chr,start,end)

#20190227 prepare rna data
rna.matrix.merged<-fread("~/Desktop/NGS_postdoc/data/RNA.expr.mat/20190519.tech.rep.collapsed.RNA.matrix.tpm.tsv")
rna.matrix.merged<-merge(rna.matrix.merged,gene.anno,by="gene") %>% merge(sample.anno[,.(cellid,stage,cnv_info)],by="cellid")
rna.matrix.merged[,logTPM:=log2(tpm+1)]
setkey(rna.matrix.merged,chr,start,end)

# filter out chr x y
rna.matrix.merged[,unique(chr)]
rna.matrix.merged<-rna.matrix.merged[chr %in% paste0("chr",c(1:19))]
# filter out qc.rna failed (rna matrix is not included )
rna.matrix.merged[cellid %in% sample.anno[qc.rna==T,cellid]]

# select genes which expressed average not less than 1 cross whole stage, in this way to have same genes set used for all cells
ref.value<-rna.matrix.merged[cnv_info=="normal",mean(logTPM),.(symbol)]
ref.value[V1>=1]
rna.matrix.merged<-merge(rna.matrix.merged,ref.value,by=c("symbol"))
tmp<-rna.matrix.merged[ V1>=1]
#tmp[,expr:=logTPM/V1]
tmp$chr<-factor(tmp$chr,levels=c(paste0("chr",c(1:19)))) #,"X","Y"

#normalized to average expression level by each gene
#res<-dcast(tmp, symbol+chr+start+end~cellid,value.var = "expr") %>% setkey(chr,start,end) %>% as.data.frame()

res<-dcast(tmp, symbol+chr+start+end~cellid,value.var = "logTPM") %>% setkey(chr,start,end) %>% as.data.frame()
#dim(res)
#res[1:10,1:6]
#table(res$chr)
#### average closest 100 gene and move to show the expression levels  
all_cnv <- lapply(split(res,res$chr), function(x){
  # x=split(res,res$chr)[[1]]
  anno=x[,1:4]
  dat=x[,5:ncol(res)]
  if(nrow(dat)>100){
    cnv <- lapply(51:(nrow(dat)-50), function(i){
      this_cnv <- unlist( lapply(1:ncol(dat), function(j){
        sum(dat[(i-50):(i+50),j])/101
      }))
      return(this_cnv)
    })
    cnv=do.call(rbind,cnv)
    cnv=cbind(anno[51:(nrow(x)-50),],cnv)
    # cnv[1:4,1:8]
  }else{
    return(NULL)
  }
})

all_cnv=do.call(rbind,all_cnv) 
#head(all_cnv[1:4,1:8])
#table(all_cnv$chr)


#library(pheatmap)
D=(t(scale(all_cnv[,5:ncol(all_cnv)] ))) 

#dim(D)
colnames(D)=paste0('genes_',1:ncol(D))
rownames(D)=colnames(res[,5:ncol(res)])
#dim(res)

tmp<-melt(D) %>% as.data.table
class(tmp)
tmp[,stage:=substring(Var1,1,4)]

tmp[,ave.gene:=mean(value),.(Var2,stage)]
tmp[,final.value:=value-ave.gene]
tmp[final.value>2, final.value:=2]
tmp[final.value< -2, final.value:= -2]
tmp<-dcast(tmp,Var1~Var2,value.var = "final.value") %>% as.data.frame() %>% tibble::column_to_rownames("Var1")
tmp[1:5,1:5]
anno_row<-sample.anno[cellid %in% rownames(tmp),.(cellid,stage,embryo_id,cnv_info,sex,GCH_sites)]%>% setkey(stage,embryo_id,cellid)%>% .[,embryo_id:=NULL]%>% as.data.frame()%>% tibble::column_to_rownames("cellid")
anno_col = data.frame(
  chr= factor(all_cnv$chr,levels =  unique(all_cnv$chr))
)
rownames(anno_col) = colnames(tmp)

library(plyr)
col.gaps<-as.data.frame(table(all_cnv$chr))
col.gaps$group<-1
class(col.gaps)
#
col.gaps<- ddply(col.gaps, "group",
                 transform, sep=cumsum(Freq))


row.gaps<-as.data.frame(table(rowanno$stage))
row.gaps$group<-1
class(row.gaps)
#library(plyr)
row.gaps<- ddply(row.gaps, "group",
                 transform, sep=cumsum(Freq))

pdf(paste0("~/Desktop/NGS_postdoc/data/CNV/20191010.RNA_CNV.tpm.scale.1cell.2sub.stage.ave.pdf"),width = 15,height=15)
p1<-pheatmap::pheatmap(tmp[rownames(anno_row),],cluster_rows = F,cluster_cols = F,
                       cellwidth=0.05,cellheight = 4,
                       #main = i,
                       # col=rev(cols),
                       color = colorRampPalette(c("royalblue","white","orangered"))(100),
                       border_color = NA,
                       gaps_col=col.gaps$sep,
                       annotation_col=anno_col,
                       annotation_row = anno_row,gaps_row = row.gaps$sep,
                       annotation_colors = annocolor,
                       show_rownames=F,show_colnames=F)
print(p1)
dev.off()




