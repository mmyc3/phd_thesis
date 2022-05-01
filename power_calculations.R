#!/usr/bin/env Rscript

## This script calculates power across a range of allele frequencies and odds ratios for a given sample size and case:control ratio

library(genpwr)
library(ggplot2)

# Power calculations
pw <- genpwr.calc(calc = "power", model = "logistic", ge.interaction = NULL,
                  N=23859, Case.Rate=0.005563282, k=NULL,
                  MAF=c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5), 
                  OR=seq(1.1,4,0.1), 
                  Alpha=0.00000005,  
                  True.Model=c("Additive"), 
                  Test.Model=c("Additive"))

a <- pw[,c(3,4,9)]
names(a)[3] <- "POWER"

# Generate graph 
png("power_gwas_graph.png", type="cairo", res=300, height=6, width=6, units="in")
ggplot(data=a, aes(x=OR, y=POWER, group=MAF, color=MAF)) + 
  geom_line(aes(color=factor(MAF))) + 
  theme_classic() +
  labs(x="Odds Ratio", y="Power") +
  scale_fill_manual(values=rainbow(7)) +
  scale_x_continuous(breaks=seq(1,max(a$OR), by=0.5)) +
  scale_y_continuous(breaks=seq(0,1,0.2))
dev.off()            
