#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=TRUE)
trait=args[1]

if (length(args)==0) {
  stop("Please supply trait.n", call.=FALSE)
}

library(data.table)
library(tidyr)
library(qqman)
library(GWASTools)

#Genomic inflation factor function
lambda<-function(pvalues){
chisq <- qchisq(1-pvalues,1)
lambda=median(chisq,na.rm=TRUE)/qchisq(0.5,1)
return(lambda)
}

#Read GWAS file
gwas<-fread(paste0("step2_v0.42.1_",trait,".SAIGE.gwas_concat_misHWEfiltered.txt"))
gwas <- unite(gwas, "SNP", c("CHR", "POS", "Allele1", "Allele2"),sep=":", remove=FALSE)
p <- length(unique(gwas$SNP))
p
#snpsOfInterest <- read.table("snps_of_interest.txt", header=F)
#snpsOfInterest <- as.vector(snpsOfInterest$V1)

# Adjust p values by FDR and Bonferroni
gwas$fdr <- p.adjust(gwas$p.value, method="fdr")
gwas$bonferroni <- p.adjust(gwas$p.value, method="bonferroni")

fdr <- gwas[gwas$fdr<0.05,]
bon <- gwas[gwas$bonferroni<0.05,]
sig <- gwas[gwas$P < 0.00000005,]
all_sig <- unique(rbind(fdr, bon, sig))
write.table(all_sig, file=paste0(trait,"_top_snps.txt"), col.names=T, row.names=F, quote=F, sep="\t")

#Change chromosome X to chr23 to be compatible with qqman
gwas$CHR=gsub("chrX","chr23",gwas$CHR)

#change chromosome variable to numeric by substituting the "chr" character prefix
gwas$CHR=as.numeric(gsub("chr","",gwas$CHR))

#Change colnames of POS and p.value to BP and P to be compatible with qqman
colnames(gwas)[which(colnames(gwas)=="POS")]="BP"
colnames(gwas)[which(colnames(gwas)=="p.value")]="P"

gwas <- gwas[!is.na(gwas$P), ]
gwas$BP <- as.numeric(gwas$BP)
gwas$P <- as.numeric(gwas$P)

#Create manhattan plot
png(paste0("saige_",trait,"_manhattan_maf0.001.png"), type="cairo", units="in", res=300, width=12, height=6)
par(omi=c(0,0.5,0,0))
manhattan(gwas, cex = 0.6, cex.axis = 0.6, col = c("lightskyblue3", "deepskyblue4"), genomewideline =-log10(0.00000005), suggestiveline = FALSE, chrlabs=c(1:22,"X"))
dev.off()

#Create qqplot with calculated inflation factor in the title
png(paste0("saige_",trait,"_qqplot.png"), type = "cairo", units="in", res=300, width=6, height=6)
par(omi=c(0,0.5,0,0))
qqPlot(gwas$P,main=paste0("lambda=",round(lambda(gwas$P),4)), cex=0.6)
dev.off()
