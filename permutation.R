#!/usr/bin/env Rscript

# Permutation testing 

args <- commandArgs(trailingOnly = TRUE)
trait <- args[1]

library(coin)
library(data.table)
library(dplyr)

bed="encode"
n=26018
svlen=0.05
pheno="puv"

ac = function(x){as.character(x)}
an = function(x){as.numeric(as.character(x))}
cc = function(x){length(unique(x[,5]))}

# Read in SV files and ped file - merge together  
vars <- c("CNV", "DEL", "DUP", "INV")
for(sv in vars){
  both<- read.delim(paste0("Rare_StructuralVariants_",bed,"_",sv,"_",n,"x_participants.txt"), header=T, sep="\t")
  both <- both[!(both$CHROM=="chrX" | both$CHROM=="chrY"),]
  both <- both[,c(8,9,10,4,1,2,3,7,11,5)]
  both$SVlengthkb <-an(both$LENGTH.Mb.)*1000
  both <- both[,-c(9)]
  both <- na.omit(both)
  both <- both[!duplicated(both[,c(1,4,5)]),]
  assign(paste0(sv,"_gencode"),both)
}

df <- rbind(CNV_gencode, DEL_gencode, DUP_gencode, INV_gencode)
both <- df

ped <- fread(paste0("~/re_gecip/renal/mchan/CAKUT/pheno/",trait,"_ancestry_matched_controls_unrelated_aggv2_pcs.ped"), data.table=F)
ped <- na.omit(ped)
cases <- as.vector(unique(ped$FID[which(ped$PHENO==1)]))
controls <- as.vector(unique(ped$FID[which(ped$PHENO==0)]))
nc <- an(length(cases))
ncon <- an(length(controls))
a <- both[which(both$Part_ID %in% cases),]
b <- both[which(both$Part_ID %in% controls),]

all_cre=c("dELS", "pELS", "PLS", "CTCF-only", "DNaseH3K4me3")

# Analyse by CRE
for (cre in all_cre){
  l <- both[grep(paste0(cre),as.vector(both$GeneSymbol)),]
  m <- merge(ped, l, by.x="FID", by.y="Part_ID", all.x=T)
  
  m$cnv <- ifelse(m$CONSEQUENCE=="CNV",1,0)
  m$inv <- ifelse(m$CONSEQUENCE=="INV",1,0)
  m$del <- ifelse(m$CONSEQUENCE=="DEL",1,0)
  m$dup <- ifelse(m$CONSEQUENCE=="DUP",1,0)
  m$sv <- ifelse(m$CONSEQUENCE=="SV",1,0)
  df <- m[,c(1,6,26:30)]

  # Reproducibility 
  set.seed(2021)
  new_ps=c()
  ps=c()

  # Loop through 5 SV types 
  for(i in 1:5){
    d <- df[,c(1,2,i)]
    names(d)[3] <- "sv"
    d <- d %>% arrange(FID, -sv) %>% filter(duplicated(FID)==FALSE)
    d$sv[is.na(d$sv)] <- 0

    # number of observations to sample
    n <- length(d$PHENO)
    # number of permutations
    P <- 10000
    # variable to resample from 
    variable <- d$PHENO

    # Calculate baseline fisher
    row_fisher <- function(m, alt = 'two.sided', cnf = 0.95) {
    f <- fisher.test(m, alternative = alt, conf.level = cnf)
    return(p_val = f$p.value)
  }

  m <- matrix(c(tabulate(d$sv[d$PHENO==1]), nc-tabulate(d$sv[d$PHENO==1]), tabulate(d$sv[d$PHENO==0]), ncon-tabulate(d$sv[d$PHENO==0])), nrow=2)
  p <- row_fisher(m)

  # initialize matrix to store permutation data
  PermSamples <- matrix(0,nrow=n, ncol=P)

  for(i in 1:P){
    PermSamples[,i] <- sample(variable, size=n, replace=FALSE)
  }

  # initialize vectors to store Test-stats
  Perm.test.stat <- rep(0,P)
  count=0

  for(i in 1:P){
    m <- matrix(c(tabulate(d$sv[PermSamples[,i]=="1"]),nc-tabulate(d$sv[PermSamples[,i]=="1"]),tabulate(d$sv[PermSamples[,i]=="0"]), ncon-tabulate(d$sv[PermSamples[,i]=="0"])), nrow=2)
    Perm.test.stat[i] <- row_fisher(m)
    count=count+1
    print(count)
  }

# Calculate P value
# ci <- quantile(Perm.test.stat, probs=c(0.025, 0.975))[1]
  new_p <- as.numeric(length(which(Perm.test.stat < p)))/P
  new_ps <- append(new_ps, new_p)
}

# Density plot of all Permutation test-stats
  plot(density(Perm.test.stat))
  abline(v=p, col="blue", lty="dotted")
  abline(v=new_p, col="red", lty="dotted")

# Using coin
# independence_test(d$PHENO ~ d$sv)
