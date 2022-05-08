#!/usr/bin/env Rscript

## This script divides the rare SV calls into types and per desired gene if needed in both the cases and the controls
# outputs a table of p values for further analysis
# INPUT: text file of MANTA and CANVAS calls intersecting with specified regions, PASS filter and sample ID.

args <- commandArgs(trailingOnly = TRUE)
pheno <- args[1]

library(ggplot2)
library(data.table)
library(dplyr)

bed="gencode"
svlen=0.05
pheno="cakut"
n=26018

ac = function(x){as.character(x)}
an = function(x){as.numeric(as.character(x))}
cc = function(x){length(unique(x[,5]))}

vars <- c("CNV", "DEL", "DUP", "INV")

for(sv in vars){
both<- read.delim(paste0("Rare_StructuralVariants_",bed,"_",sv,"_",n,"x_participants.txt"), header=T, sep="\t")
both <- both[!(both$CHROM=="chrX" | both$CHROM=="chrY"),]
both <- both[,c(8,9,10,4,1,2,3,7,11,5)]
both$SVlengthkb <-an(both$LENGTH.Mb.)*1000
both <- both[,-c(9)]
both <- na.omit(both)
both <- both[!duplicated(both[,c(1,4,5)]),]
assign(paste0(sv,"_exon"),both)
}

df <- rbind(CNV_exon, DEL_exon, DUP_exon, INV_exon)
both <- df

ped <- fread(paste0("/re_gecip/renal/mchan/CAKUT/pheno/",pheno,"_ancestry_matched_controls_unrelated_aggv2_pcs.ped"), data.table=F)
ped <- na.omit(ped)
cases <- as.vector(unique(ped$FID[which(ped$PHENO==1)]))
controls <- as.vector(unique(ped$FID[which(ped$PHENO==0)]))
nc <- an(length(cases))
ncon <- an(length(controls))
cases_filtered <- both[which(both$Part_ID %in% cases),]
controls_filtered <- both[which(both$Part_ID %in% controls),]

# Read in list of GENCODE v29 protein coding genes (n=19907)
panel <- read.table("/re_gecip/renal/mchan/CAKUT/lists/gencode_v29.bed", header=F)
panel <- as.vector(unique(panel$V1))

# Create empty vectors
del_fisher=c()
dup_fisher=c()
cnv_fisher=c()
inv_fisher=c()
sv_fisher=c()

# Initialise counter
count=0

for(gene in panel){
  a <- subset(cases_filtered, grepl(paste("\\b",gene,"\\b", sep=""), GeneSymbol))
  b <- subset(controls_filtered, grepl(paste("\\b",gene,"\\b", sep=""), GeneSymbol))

  count=count+1
  print(count)

#subset into SV type
del <- subset(a, a$CONSEQUENCE == "DEL" & a$SVlengthkb >= svlen)
dup <- subset(a, a$CONSEQUENCE == "DUP" & a$SVlengthkb >= svlen)
cnv <- subset(a, a$CONSEQUENCE == "CNV" & a$SVlengthkb >= svlen)
inv <- subset(a, a$CONSEQUENCE == "INV" & a$SVlengthkb >= svlen)
all <- rbind(cnv,del,dup,inv)

del_b <- subset(b, b$CONSEQUENCE == "DEL" & b$SVlengthkb >= svlen)
dup_b <- subset(b, b$CONSEQUENCE == "DUP" & b$SVlengthkb >= svlen)
cnv_b <- subset(b, b$CONSEQUENCE == "CNV" & b$SVlengthkb >= svlen)
inv_b <- subset(b, b$CONSEQUENCE == "INV" & b$SVlengthkb >= svlen)
all_b <- rbind(cnv_b,del_b, dup_b,inv_b)

##create data frame for no individuals with rare (MAF< 0.001) SV >= 50bp and run fishers exact test on each row
DEL <- c(an(length(unique(del$Part_ID))),nc-an(length(unique(del$Part_ID))),an(length(unique(del_b$Part_ID))), ncon-an(length(unique(del_b$Part_ID))))
DUP <- c(an(length(unique(dup$Part_ID))),nc-an(length(unique(dup$Part_ID))),an(length(unique(dup_b$Part_ID))), ncon-an(length(unique(dup_b$Part_ID))))
CNV <- c(an(length(unique(cnv$Part_ID))),nc-an(length(unique(cnv$Part_ID))),an(length(unique(cnv_b$Part_ID))), ncon-an(length(unique(cnv_b$Part_ID))))
INV <- c(an(length(unique(inv$Part_ID))),nc-an(length(unique(inv$Part_ID))),an(length(unique(inv_b$Part_ID))), ncon-an(length(unique(inv_b$Part_ID))))
TOTAL <- c(an(length(unique(all$Part_ID))),nc-an(length(unique(all$Part_ID))),an(length(unique(all_b$Part_ID))),ncon-an(length(unique(all_b$Part_ID))))
d <-t(data.frame(CNV, DEL, DUP, INV, TOTAL))
colnames(d) <- c("cases_aff", "cases_unaff", "controls_aff", "controls_unaff")
d <- as.data.frame(d)

## Set Fishers exact function
  row_fisher <- function(row, alt = 'two.sided', cnf = 0.95) {
    f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
    return(c(row,
             p_val = f$p.value,
             or = f$estimate[[1]],
             or_ll = f$conf.int[1],
             or_ul = f$conf.int[2]))
  }
  d <- data.frame(t(apply(d, 1, row_fisher)))

  del_fisher <- append(del_fisher, d[2,5])
  dup_fisher <- append(dup_fisher, d[3,5])
  cnv_fisher <- append(cnv_fisher, d[1,5])
  inv_fisher <- append(inv_fisher, d[4,5])
  sv_fisher <- append(sv_fisher, d[5,5])
}

df <- data.frame(panel, cnv_fisher, del_fisher, dup_fisher, inv_fisher, sv_fisher)
rownames(df) <- df[,1]
colnames(df) <- c("GENE", "CNV", "DEL", "DUP", "INV", "ALL_SV")
write.table(df, paste0(trait,"_exome_wide_SV_summary_stats.txt"), col.names=T, row.names=F, sep="\t", quote=F)
