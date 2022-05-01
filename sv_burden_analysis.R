#!/usr/bin/env Rscript

## This script divides rare (MAF < 0.1%) SV calls into types (CNV, INV, DEL, DUP as per MANTA/CANVAS) and per cCRE (dELS, pELS, PLS, CTCF-only, DNase-H3K4me3)
## Performs a Fisher's exact test to compare the burden between cases and controls
## Inputfile is a tab-delimited text file listing all rare (MAF < 0.1%) SVs (separated by type) within the cohort of interest annotated with the genes or cCREs they intersect.
## Phenofile is a tab-delimited text file in the format FID, IID, PATID, MATID, SEX, PHENO (1=affected, 0-unaffected)

library(ggplot2)
library(data.table)
library(dplyr)
library(Cairo)

## Choose SV length threshold (in kb)
svlen=0.05
input="inputfile"
pheno="phenofile"

ac = function(x){as.character(x)}
an = function(x){as.numeric(as.character(x))}
cc = function(x){length(unique(x[,5]))}

## Read in input file per SV type, remove chrX and chrY and duplicated SVs

vars <- c("CNV", "DEL", "DUP", "INV")

for(sv in vars){
  both <- fread(paste0(sv,"_",input))
  both <- both[!(both$CHROM=="chrX" | both$CHROM=="chrY"),]
  both <- both[,c(8,9,10,4,1,2,3,7,11,5)]
  both$SVlengthkb <-an(both$LENGTH.Mb.)*1000
  both <- both[,-c(9)]
  both <- na.omit(both)
  both <- both[!duplicated(both[,c(1,4,5)]),]
  assign(paste0(sv),both)
}

## Combine in dataframe and extract those SVs > svlen
y <- rbind(CNV, DEL, DUP, INV)
y <- y[y$SVlengthkb > svlen,]

## Read in phenotype file and separate SVs seen in cases and controls
ped <- fread(paste0(pheno), data.table=F)
ped <- na.omit(ped)
cases <- as.vector(unique(ped$FID[which(ped$PHENO==1)]))
controls <- as.vector(unique(ped$FID[which(ped$PHENO==0)]))
nc <- an(length(cases))
ncon <- an(length(controls))
a <- df[which(df$Part_ID %in% cases),]
b <- df[which(df$Part_ID %in% controls),]

# Perform Fisher's exact per SV type

svlist=list()

for(sv in vars){
  SV <- subset(a, a$CONSEQUENCE == sv & a$SVlengthkb >= svlen)
  SV_b <- subset(b, b$CONSEQUENCE == sv & b$SVlengthkb >= svlen)

  #create data frame for no individuals with rare (MAF< 0.001) SV >= 50bp and run fishers exact test on each row
  SV_c <- c(an(length(unique(SV$Part_ID))),nc-an(length(unique(SV$Part_ID))),an(length(unique(SV_b$Part_ID))), ncon-an(length(unique(SV_b$Part_ID))))
  d <-t(data.frame(SV_c))
  colnames(d) <- c("cases_aff", "cases_unaff", "controls_aff", "controls_unaff")
  d <- as.data.frame(d)

  row_fisher <- function(x) {
    mat <- matrix(as.numeric(x[1:4]), nrow = 2)
    f <- fisher.test(mat, alternative = "two.sided", conf.level = 0.95)
    return(c(x,
             p_val = f$p.value,
             or = f$estimate[[1]],
             or_ll = f$conf.int[1],
             or_ul = f$conf.int[2]))
  }

  ci <- function(x) {
    mat <- matrix(as.numeric(x[1],x[3],nc,ncon), nrow=2)
    c <- prop.test(x[1],nc, alternative="two.sided", conf.level=0.95)
    h <- prop.test(x[3],ncon, alternative="two.sided", conf.level=0.95)
    return(c(x,ci_up=c$conf.int[2]*100,
             ci_low=c$conf.int[1]*100,
             ci_up_con=h$conf.int[2]*100,
             ci_low_con=h$conf.int[1]*100))
  }

  d <- as.data.frame(t(apply(d, 1, row_fisher)))
  d <- as.data.frame(t(apply(d, 1, ci)))
  d$prop <- an(d[,1]/nc*100)
  d$prop_con <- an(d[,3]/ncon*100)
  d$SV <- sv
  svlist[[sv]] <- d
}

df <- data.table::rbindlist(svlist)

# Compare size of SVs between cases and controls
for(sv in vars){
  SV <- subset(a, a$CONSEQUENCE == sv & a$SVlengthkb >= svlen)
  SV_b <- subset(b, b$CONSEQUENCE == sv & b$SVlengthkb >= svlen)
  SV$cohort <- "CASE"
  SV_b$cohort <- "CONTROL"
  SV_all <- rbind(SV[,c(11,10)], SV_b[,c(11,10)])
  print(paste("Cases median is ",median(SV$SVlengthkb),"kb",sep=""))
  print(paste("Cases IQR is ",IQR(SV$SVlengthkb),"kb",sep=""))
  print(paste("Controls median is ",median(SV_b$SVlengthkb),"kb",sep=""))
  print(paste("Controls IQR is ",IQR(SV_b$SVlengthkb),"kb",sep=""))
  print(try(wilcox.test(SVlengthkb~cohort, data=SV_all, exact=FALSE)))
}

# Analyze by cCRE type

all_cre <- c("dELS", "PLS", "pELS", "CTCF-only", "DNase-H3K4me3")

datalist=list()
by_cre=list()

for (cre in all_cre){
  # Analyse by CRE
  s <- a[grep(paste0(cre),as.vector(a$GeneSymbol)),]
  t <- b[grep(paste0(cre),as.vector(b$GeneSymbol)),]

  for(sv in vars){
    SV <- subset(s, s$CONSEQUENCE == sv & s$SVlengthkb >= svlen)
    SV_b <- subset(t, t$CONSEQUENCE == sv & t$SVlengthkb >= svlen)

    #create data frame for no individuals with rare (MAF< 0.001) SV >= 50bp and run fishers exact test on each row
    SV_c <- c(an(length(unique(SV$Part_ID))),nc-an(length(unique(SV$Part_ID))),an(length(unique(SV_b$Part_ID))), ncon-an(length(unique(SV_b$Part_ID))))
    d <-t(data.frame(SV_c))
    colnames(d) <- c("cases_aff", "cases_unaff", "controls_aff", "controls_unaff")
    d <- as.data.frame(d)
    d <- as.data.frame(t(apply(d, 1, row_fisher)))
    d <- as.data.frame(t(apply(d, 1, ci)))
    d$prop <- an(d[,1]/nc*100)
    d$prop_con <- an(d[,3]/ncon*100)
    d$cre <- cre
    d$SV <- sv
    datalist[[sv]] <- d
  }
  
  by_cre[[cre]] <- data.table::rbindlist(datalist)
}

all <- data.table::rbindlist(by_cre)

# Create graph for CRE
cohort <- rep(c("CONTROL", "PUV"),20)
freq <- as.vector(t(cbind(all$prop_con,all$prop)))
ci_up <- as.vector(t(cbind(all$ci_up_con, all$ci_up)))
ci_low <- as.vector(t(cbind(all$ci_low_con, all$ci_low)))
sv <- rep(c("CNV", "CNV","DEL", "DEL", "DUP", "DUP", "INV","INV"),5)
cre <- c(rep("dELS",8), rep("pELS",8),rep("PLS",8), rep("CTCF-only",8), rep("DNase-H3K4me3",8))
q <- data.frame("COHORT"=cohort, "SV"=sv, "cCRE"=cre, "FREQ"=freq, "CI_UP"=ci_up, "CI_LOW"=ci_low)
q$cCRE <- factor(q$cCRE,levels=unique(q$cCRE))

dodge <- position_dodge(width = 0.9)
col=c("deepskyblue4", "darkorange", "darkolivegreen4", "darkorchid4")

k <- ggplot(q, aes(x = cCRE, y = FREQ, fill=SV, alpha=COHORT)) +
  geom_bar(stat = "identity", position='dodge', width=0.8) +
  scale_alpha_manual(values=c(0.5,1)) +
  scale_fill_manual(values=col) +
  geom_errorbar(aes(ymax = CI_UP, ymin = CI_LOW), position = position_dodge(0.8), width = 0.1) +
  theme_classic() +
  facet_grid(.~SV) +
  theme(text=element_text(size=10), axis.text.x=element_text(angle=45, hjust=1)) +
  scale_y_continuous(breaks=seq(0,100,10), limits=c(0,100)) +
  labs(x="Candidate cis-regulatory element (cCRE)", y="Proportion of individuals\nwith \u22651 SV (%)")

CairoPDF("SV_cre_burden.pdf", width=10)
k
dev.off()
