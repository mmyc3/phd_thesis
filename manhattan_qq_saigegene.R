#!/usr/bin/env Rscript

# This script generates a Manhattan and Q-Q plot using the results from SAIGE-GENE 

library(tidyr)
library(data.table)
library(qqman)

args <- commandArgs(trailingOnly=TRUE)
trait=args[1]

# Read in SAIGE-GENE results text file
setwd(paste0("/re_gecip/renal/mchan/CAKUT/saige_gene/",trait))
a <- fread(paste0("all_chr_",trait,"_maf0.001_cadd20.SAIGE.gene.txt"), data.table=F)
a <- separate(a, col="Gene", into=c("Gene", "ENS"), sep="\\_")

# calculate FDR p values for reference
a$fdr <- p.adjust(a$Pvalue, method="fdr")
sig <- a[a$fdr<0.05,]
b <- a[,c(1,3)]

# Read in list of GENCODE gene names and positions 
gene <- fread("/re_gecip/renal/mchan/CAKUT/lists/gencode_v29_genenames.txt", header=T, data.table=F)
gene <- unique(gene)
gene <- gene[!duplicated(gene$name2),]
names(gene) <- c("CHR", "START", "END", "Gene")
gene$START <- as.numeric(gene$START)
gene$END <- as.numeric(gene$END)
gene$BP <- (gene$START+gene$END)/2

# Merge results and gene file, rename columns and chrX
c <- merge(b, gene, by="Gene", all.x=T)
c <- c[,c(1:3,6)]
names(c) <- c("SNP", "P", "CHR", "BP")
c$CHR=gsub("chrX","chr23",c$CHR)
c$CHR=as.numeric(gsub("chr","",c$CHR))
c <- c[!is.na(c$P), ]
c$BP <- as.numeric(c$BP)
c$P <- as.numeric(c$P)
d <- na.omit(c)

# How many genes tested 
p=length(unique(d$SNP))

# Plot Manhattan
png(paste0(trait,"_saigegene_manhattan_maf0.001.png"), type="cairo", width=12, height=6, units="in", res=300)
par(omi=c(0,0.5,0,0))
manhattan(d, ylim=c(0,8), cex.axis=0.6, col = c("lightskyblue3", "deepskyblue4"), annotatePval=0.0001, genomewideline =-log10(0.05/length(unique(d$SNP))), suggestiveline = FALSE, chrlabs=c(1:22, "X"))
dev.off()

# Create qqplot
png(paste0(trait,"saigegene_qqplot.png"), type = "cairo", units="in", res=300, width=6, height=6)
par(omi=c(0,0.5,0,0))
qqPlot(gwas$P)
dev.off()
