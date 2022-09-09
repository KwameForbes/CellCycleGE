library(DESeq2)
library(tximport)
library(dplyr)
install.packages('scrime')
library(DESeq2)
library("scrime")


#This is the pathway to the files that undergone salmon quantification
pathway <- "/pine/scr/k/w/kwamek/pengda_collab/Quantfiles2"

dir <- list.files(path = pathway, pattern = NULL, all.files = FALSE,
                  full.names = TRUE, recursive = FALSE,
                  ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

#list of files in quants.sf files
sample.quants <- file.path(dir,"quant.sf")
##building txi
t2g <- read.table(file = "/proj/RNA_lab/ZARD/RIP/tx2genes.txt", sep = "\t", header = T)
head(t2g)
t2g <- t2g[,c("ENST","ENSG","Gene_name")]
dim(t2g)
genes <- t2g[,c("ENSG","Gene_name")]
genes <- unique(genes, by="ENSG")
head(genes)



txi <- tximport(sample.quants,type="salmon",tx2gene = t2g,lengthCol = "Length",existenceOptional=T,ignoreAfterBar = T)

sample.data <- read.table("/pine/scr/k/w/kwamek/pengda_collab/TimePoints.txt",header = F)
name <- sample.data$V1
TimePoints <- sample.data$V2
samples <- cbind(name,TimePoints)

#Building deseq object
dds <- DESeqDataSetFromTximport(txi,colData = samples,design = ~ TimePoints)
dds <- DESeq(dds)

#library size normalization
GECellCycle = as.data.frame(counts(dds,normalize =TRUE))
GECellCycle = rename_with(GECellCycle, ~ (gsub("V", "TP_", .x,fixed = TRUE)))
#write.csv(GECellCycle,"GECellCycle.csv")
#library size normalized and row wise normalized
GECellCyclescaled = t(apply(GECellCycle, 1, function(x)(x-min(x))/(max(x))))
GECellCyclescaled = as.data.frame(GECellCyclescaled)
GECellCyclescaled = na.omit(GECellCyclescaled)
#write.csv(GECellCyclescaled,"GECellCyclescaled.csv")

#pheatmap library size normalized z scored
pheatmap(GECellCycle, show_rownames=FALSE,cluster_cols= F,cluster_rows = F,scale = "row",
         cellwidth = 10,color=colorRampPalette(c("red", "white", "blue"))(150),
         border_color = NA)

