#!/usr/bin/env Rscript

library(Rlabkey)
library(dplyr)
library(ggplot2)
library(Cairo)

labkey.setDefaults(baseUrl="https://labkey-embassy.gel.zone/labkey/")

setwd("~/re_gecip/renal/mchan/CAKUT/clinical")

# Pull out all CAKUT probands with HPO terms
query <- "SELECT DISTINCT rd.participant_id, rd.plate_key, rd.rare_diseases_family_id, rd.biological_relationship_to_proband, rd.genome, rd.path, rd.inferred_sex_karyotype, hp.normalised_hpo_id, hp.normalised_hpo_term, rd.genome_build, par.year_of_birth, par.participant_phenotypic_sex, par.consanguinity, par.father_affected, par.mother_affected, par.full_brothers_affected, par.full_sisters_affected, par.participant_ethnic_category
          FROM rare_disease_analysis AS rd
          LEFT JOIN genome_file_paths_and_types AS sq
          ON rd.participant_id = sq.participant_id
          LEFT JOIN rare_diseases_participant_phenotype AS hp
          ON rd.participant_id = hp.participant_id
          LEFT JOIN participant AS par
          ON rd.participant_id = par.participant_id
          WHERE rd.normalised_specific_disease = 'Congenital Anomaly of the Kidneys and Urinary Tract (CAKUT)'
          AND hp.hpo_present = 'Yes'"

mysql <- labkey.executeSql(
  schemaName="lists",                                                 # Do not change this
  colNameOpt = "rname",                                               # Do not change this
  maxRows = 100000000,                                                # Do not change this
  folderPath="/main-programme/main-programme_v10_2020-09-03",          # This can be changed to different main programme releases
  sql = query                                                         # This can be changed to your query of choice
)

x <-data.frame(mysql)
y <- unique(x)
y <- y[!(y$plate_key=="N/A"),]

# Add columns for age, sex, family history, consanguinity and extra-renal features 
extra <- read.delim("extra_renal_hpo_list", header=F)
y$extra_renal <- ifelse(y$normalised_hpo_term %in% extra$V1, 1, 0)
y$SEX <- ifelse(y$participant_phenotypic_sex=="Male",1,2)
y$age <- 2021-y$year_of_birth
age <- unique(y[,c(1,15)])
summary(age$age)

y %>% group_by(participant_phenotypic_sex) %>% summarise(count=n_distinct(participant_id))
w <- y %>% group_by(participant_id) %>% filter(extra_renal == 1)
v <- y %>% group_by(participant_id) %>% filter(extra_renal == 0)
t <- w$participant_id
v <- v[!(v$participant_id %in% t),]
u <- rbind(w,v)
u %>% group_by(extra_renal) %>% summarise(count=n_distinct(participant_id))
y %>% group_by(participant_id) %>% summarise(count=n_distinct(normalised_hpo_term))
y[is.na(y)] <- 0
y$fh <- ifelse((y$father_affected=="Yes" | y$mother_affected=="Yes" | y$full_brothers_affected!=0 | y$full_sisters_affected!=0),1,0)
y %>% group_by(fh) %>% summarise(count=n_distinct(plate_key))
y$cons <- ifelse((y$consanguinity=="Yes" | y$consanguinity=="Possible"),1,0)
y %>% group_by(cons) %>% summarise(count=n_distinct(plate_key))

# ethnicity breakdown 
y %>% filter(grepl("^White", participant_ethnic_category)) %>% summarise(count=n_distinct(plate_key))
y %>% filter(grepl("^Asian", participant_ethnic_category)) %>% summarise(count=n_distinct(plate_key))
y %>% filter(grepl("^Black", participant_ethnic_category)) %>% summarise(count=n_distinct(plate_key))
y %>% filter(grepl("Chinese", participant_ethnic_category)) %>% summarise(count=n_distinct(plate_key))
y %>% filter(grepl("^Mixed", participant_ethnic_category)) %>% summarise(count=n_distinct(plate_key))
y %>% filter(grepl("^Other", participant_ethnic_category)) %>% summarise(count=n_distinct(plate_key))
y %>% filter(grepl("^Not", participant_ethnic_category)) %>% summarise(count=n_distinct(plate_key))

# calculate numbers for each phenotype
ka <- readLines("~/re_gecip/renal/mchan/CAKUT/lists/ka_hpo.txt")
ou <- readLines("~/re_gecip/renal/mchan/CAKUT/lists/ou_hpo.txt")
cystic <- readLines("~/re_gecip/renal/mchan/CAKUT/lists/cystic_hpo.txt")
ek <- readLines("~/re_gecip/renal/mchan/CAKUT/lists/ectopic_hpo.txt")
dup <- readLines("~/re_gecip/renal/mchan/CAKUT/lists/duplex_hpo.txt")

y %>% filter(grepl(paste(ka,collapse="|"), normalised_hpo_term, ignore.case=TRUE)) %>% summarise(count=n_distinct(participant_id))
y %>% filter(grepl(paste(ou,collapse="|"), normalised_hpo_term, ignore.case=TRUE)) %>% summarise(count=n_distinct(participant_id))
y %>% filter(grepl(paste(cystic,collapse="|"), normalised_hpo_term, ignore.case=TRUE)) %>% summarise(count=n_distinct(participant_id))
y %>% filter(grepl(paste(ek,collapse="|"), normalised_hpo_term, ignore.case=TRUE)) %>% summarise(count=n_distinct(participant_id))
y %>% filter(grepl(paste(dup,collapse="|"), normalised_hpo_term, ignore.case=TRUE)) %>% summarise(count=n_distinct(participant_id))
y %>% filter(grepl("posterior urethral valve", normalised_hpo_term, ignore.case=TRUE)) %>% summarise(count=n_distinct(participant_id))
y %>% filter(grepl("vesicoureteral reflux", normalised_hpo_term, ignore.case=TRUE)) %>% summarise(count=n_distinct(participant_id))
y %>% filter(grepl("bladder exstrophy", normalised_hpo_term, ignore.case=TRUE)) %>% summarise(count=n_distinct(participant_id))

# Count number with extra-renal phenotypes 
neuro <- c("neurodevelopmental", "developmental delay", "learning disability", "intellectual disability", "motor delay", "delayed social development", "delayed speech", "autistic behavior", "autism")
y %>% filter(grepl(paste(neuro,collapse="|"), normalised_hpo_term, ignore.case=TRUE)) %>% summarise(count=n_distinct(plate_key))
cardiac <- c("abnormal heart morphology", "congenital malformation of the left heart", "cardiac")
y %>% filter(grepl(paste(cardiac,collapse="|"), normalised_hpo_term, ignore.case=TRUE)) %>% summarise(count=n_distinct(plate_key))
y %>% filter(grepl("hearing", normalised_hpo_term, ignore.case=TRUE)) %>% summarise(count=n_distinct(plate_key))
y %>% filter(grepl("gout|hyperuricemia", normalised_hpo_term, ignore.case=TRUE)) %>% summarise(count=n_distinct(plate_key))

# ESKD
hes <- final_hes_combined[final_hes_combined$participant_id %in% y$participant_id,]
hes <- unique(hes)
length(unique(hes$participant_id))
summary(hes$age_of_esrf)

# Pull out GMC exit questionnaire
query <- "SELECT rd.participant_id, rd.plate_key, par.year_of_birth, par.participant_ethnic_category, gc.case_solved_family, gc.assembly, gc.chromosome, gc.position, gc.reference, gc.alternate, gc.gene_name, gc.acmg_classification, gc.clinical_utility, gc.additional_comments, gc.publications, par.consanguinity, par.father_affected, par.mother_affected, par.full_brothers_affected, par.full_sisters_affected, hp.hpo_present, hp.hpo_id, hp.normalised_hpo_term
          FROM rare_disease_analysis AS rd
          LEFT JOIN rare_diseases_participant_phenotype AS hp
          ON rd.participant_id = hp.participant_id
          LEFT JOIN gmc_exit_questionnaire AS gc
          ON rd.participant_id = gc.participant_id
          LEFT JOIN participant as par
          ON rd.participant_id = par.participant_id
          WHERE rd.normalised_specific_disease = 'Congenital Anomaly of the Kidneys and Urinary Tract (CAKUT)'
          AND rd.participant_type = 'Proband'"

mysql <- labkey.executeSql(
  schemaName="lists",                                                 # Do not change this
  colNameOpt = "rname",                                               # Do not change this
  maxRows = 100000000,                                                # Do not change this
  folderPath="/main-programme/main-programme_v12_2021-05-06",          # This can be changed to different main programme releases
  sql = query                                                         # This can be changed to your query of choice
)

x <-data.frame(mysql)
gmc <- unique(x)
gmc <- gmc[complete.cases(gmc[,1]),]

# Add whether solved, age, FH, consanguinity 
gmc$fh <- ifelse((gmc$father_affected=="Yes" | gmc$mother_affected=="Yes" | gmc$full_brothers_affected!=0 | gmc$full_sisters_affected!=0),1,0)
gmc$solved <- ifelse((gmc$case_solved_family=="yes" | gmc$case_solved_family=="partially"),1,0)
gmc$cons <- ifelse((gmc$consanguinity=="Yes" | gmc$consanguinity=="Possible"),1,0)
gmc$age <- 2021-gmc$year_of_birth
a <- unique(gmc[,c(2,4,24:27)])
a[is.na(a)] <- 0

# output list of all solved CAKUT probands 
both <- merge(y,a, by="plate_key", all.y=TRUE)
write.table(both, "cakut_probands_b38_b37_solved.txt", sep="\t", col.names=T, row.names=F, quote=F)

# Look at exomiser for those without GMC exit questionnaire
both <- unique(y$plate_key[!y$plate_key %in% both$plate_key])

# Pull out exomiser variants
query <- "SELECT ex.participant_id, rd.plate_key, ex.assembly, ex.genomic_feature_hgnc, ex.chromosome, ex.position, ex.reference, ex.alternate, ex.genotype, ex.mode_of_inheritance, ex.penetrance, ex.rank, ex.variant_score, ex.consequence, ex.hgvs, ex.max_frequency
          FROM rare_disease_analysis AS rd
          LEFT JOIN exomiser AS ex
          ON rd.participant_id = ex.participant_id
          WHERE rd.normalised_specific_disease = 'Congenital Anomaly of the Kidneys and Urinary Tract (CAKUT)'
          AND rd.participant_type = 'Proband'
          AND ex.rank <= 5"

mysql <- labkey.executeSql(
  schemaName="lists",                                                 # Do not change this
  colNameOpt = "rname",                                               # Do not change this
  maxRows = 100000000,                                                # Do not change this
  folderPath="/main-programme/main-programme_v12_2021-05-06",          # This can be changed to different main programme releases
  sql = query                                                         # This can be changed to your query of choice
)

x <-data.frame(mysql)
exomiser <- unique(x)

# Search known CAKUT genes for exomiser prioritized variants 
genes <- readLines("../lists/cakut_genes.txt")
ex <- exomiser %>% filter(grepl(paste(genes,collapse="|"),genomic_feature_hgnc))
k <- not_ex_hpo %>% select(participant_id.x, plate_key, participant_phenotypic_sex, participant_ethnic_category.x, age, fh, solved, cons, genomic_feature_hgnc, chromosome, position, reference, alternate, genotype, mode_of_inheritance, consequence, hgvs, assembly, normalised_hpo_term) %>%
  group_by(plate_key) %>% 
  summarise_all(~toString(unique(.)))

# Calculate proportion cases solved 
as.numeric(length(which(both$solved==1)))/as.numeric(length(unique(both$participant_id)))*100
as.numeric(length(which(both$fh==1)))/as.numeric(length(unique(both$plate_key)))*100
as.numeric(length(which(both$extra_renal==1)))/as.numeric(length(unique(both$plate_key)))*100
table(both$SEX)
summary(both$age)
table(both$cons)
both$SEX <- ifelse(both$SEX==2,0,1)

# Use logisitic regression to determine predictors of a genetic diagnosis  
a <- glm(data=both, formula=solved ~ fh + cons + age + extra_renal + SEX, family="binomial")
sink("cakut_clinical_logregression")
summary(a)
or <- exp(cbind(OR=coef(a), confint(a)))
round(or,digits=1)
sink()

# HPO term breakdown
w <- y %>% group_by(normalised_hpo_term) %>% summarise(Freq=n_distinct(participant_id))
w <- w[order(-w$Freq),]
write.table(w, "cakut_hpo_breakdown_v10.txt", col.names=T, row.names=F, quote=F, sep="\t")

w$Freq <- as.numeric(w$Freq)
v <- w[1:20,]
v$normalised_hpo_term <- factor(v$normalised_hpo_term, levels=v$normalised_hpo_term[order(v$Freq)])
p <- ggplot(v, aes(x=normalised_hpo_term, y=Freq)) +
  geom_bar(stat="identity", fill="deepskyblue4") +
  coord_flip() +
  geom_text(aes(label=Freq), hjust=1.5, color="white", size=3) +
  theme_minimal() +
  labs(y="Frequency", x="HPO term")

CairoPDF("cakut_top20_hpo.pdf", width=10)
p
dev.off()

# Compare solved cases by age, FH, consanguininty, extra-renal features 
solved <- both[both$solved==1,]
nc <- as.numeric(length(unique(solved$plate_key)))
unsolved <- both[both$solved!=1,]
ncon <- as.numeric(length(unique(unsolved$plate_key)))
ethnic <- c("White: British", "White: Any other White background", "White: Irish")

sex <- c(as.numeric(length(unique(solved$plate_key[which(solved$SEX==1)]))), as.numeric(length(unique(unsolved$plate_key[which(unsolved$SEX==1)]))))
fh <- c(as.numeric(length(unique(solved$plate_key[which(solved$fh==1)]))),as.numeric(length(unique(unsolved$plate_key[which(unsolved$fh==1)]))))
er <- c(as.numeric(length(unique(solved$plate_key[which(solved$extra_renal==1)]))),as.numeric(length(unique(unsolved$plate_key[which(unsolved$extra_renal==1)]))))
cons <- c(as.numeric(length(unique(solved$plate_key[which(solved$cons==1)]))),as.numeric(length(unique(unsolved$plate_key[which(unsolved$cons==1)]))))
eth <- c(as.numeric(length(unique(solved$plate_key[which(solved$participant_ethnic_category.x %in% ethnic)]))),as.numeric(length(unique(unsolved$plate_key[which(unsolved$participant_ethnic_category.x %in% ethnic)]))))
df <- as.data.frame(rbind(sex, fh, er, cons, eth))
colnames(df) <- c("solved", "unsolved")
df$a <- nc-df$solved
df$b <- ncon-df$unsolved
df <- df[,c(1,3,2,4)]

row_fisher <- function(row, alt = 'two.sided', cnf = 0.95) {
  f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
  return(c(row,
           p_val = f$p.value,
           or = f$estimate[[1]],
           or_ll = f$conf.int[1],
           or_ul = f$conf.int[2]))
}
p <- data.frame(t(apply(df, 1, row_fisher)))

p$solved_prop <- round(p$solved/nc*100, digits=2)
p$unsolved_prop <- round(p$unsolved/ncon*100,digits=2)

sink("solved_fisher")
p
summary(solved$age, na.rm=TRUE)
summary(unsolved$age, na.rm=TRUE)
wilcox.test(age~solved, data=b, conf.int=TRUE, conf.level=0.95, exact=FALSE, na.action=na.omit)
sink()

# Different genes identified
gene <- gmc[(gmc$case_solved_family=="yes"|gmc$case_solved_family=="partially"),c(2,11,12)]
gene<-unique(gene)
gene <- na.omit(gene)
Cystic <- c("PKD1", "PKD2", "PKHD1", "DNAJB11")
Syndromic <- c("MID1", "SHANK3", "WDR35", "USH2A", "ACTB", "SOX5", "PUF60", "TRRAP", "PTPN11", "SETD2", "SF3B4")
COL4A <- c("COL4A3", "COL4A4", "COL4A5")
Ciliopathy <- c("NPHP3", "CEP290")
Other <- c("RMND1", "CLCN5")
CAKUT <- c("PBX1", "GREB1L", "TBX6", "CHD7", "PAX2")
Tubular <- c("SLC3A1", "CLCN5")

n <- c()
for(g in type){
a <- subset(gene, grepl(paste("\\b",g,"\\b", sep="", collapse="|"), gene_name))
m <- as.numeric(length(unique(a$plate_key)))
n <- append(n,m)
}

# Pie chart of genes 
type <- c("Cystic", "Syndrome", "COL4A", "Ciliopathy", "Tubular", "CAKUT", "Other")
p <- data.frame(cbind(cohort=type,freq=n))
p$freq <- as.numeric(p$freq)

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    axis.text=element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

library(scales)
library(Cairo)
pie <- ggplot(p,aes(x="",y=freq,fill=cohort)) + 
  geom_col(width=1, position='stack') +
  coord_polar("y")+
  scale_fill_brewer(palette="Dark2") +
  theme(axis.text.x=element_blank()) +
  blank_theme +
  labs(fill="Genes") +
  geom_text(aes(label=paste0(freq,"\n(",scales::percent(freq/sum(freq)),")"),x=1.2),
            position=position_stack(vjust=0.5))


CairoPDF("cakut_solved_genes_pie.pdf")
pie
dev.off()

# Diagnostic yield by phenotype
ka <- read.table("~/re_gecip/renal/mchan/CAKUT/lists/ka.txt")
table(a$solved[which(a$plate_key %in% ka$V1)])
ou <- read.table("~/re_gecip/renal/mchan/CAKUT/lists/ou.txt")
table(a$solved[which(a$plate_key %in% ou$V1)])
vur <- read.table("~/re_gecip/renal/mchan/CAKUT/lists/vur.txt")
table(a$solved[which(a$plate_key %in% vur$V1)])
cystic <- read.table("~/re_gecip/renal/mchan/CAKUT/lists/cystic.txt")
table(a$solved[which(a$plate_key %in% cystic$V1)])
puv <- read.table("~/re_gecip/renal/mchan/CAKUT/lists/puv.txt")
table(a$solved[which(a$plate_key %in% puv$V1)])
ek <- read.table("~/re_gecip/renal/mchan/CAKUT/lists/ek_hk.txt")
table(a$solved[which(a$plate_key %in% ek$V1)])
be <- read.table("~/re_gecip/renal/mchan/CAKUT/lists/bladderexstrophy.txt")
table(a$solved[which(a$plate_key %in% be$V1)])
dup <- read.table("~/re_gecip/renal/mchan/CAKUT/lists/duplex.txt")
table(a$solved[which(a$plate_key %in% dup$V1)])

library(viridis)
col=rainbow(8)
phe <- c("Kidney anomaly", "Obstructive uropathy", "Cystic dysplasia", "VUR", "Ectopic kidney", "Duplex kidney", "PUV", "Bladder exstrophy")
per <-c(5.9,3.9,10.7,0.9,2.8,0,0,0)
t <- data.frame(Phenotype=phe, Proportion=per)
t$Phenotype <- factor(t$Phenotype, levels = t$Phenotype)
p <-ggplot(t, aes(x=Phenotype, y=Proportion, fill=Phenotype))+
  geom_bar(stat="identity") +
  theme_classic() +
  geom_text(aes(label=Proportion), vjust=1.5, color="white", size=4) +
  scale_fill_brewer(palette="Dark2") +
  labs(x="Phenotype", y="Diagnostic yield (%)") +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  theme(legend.position="none") +
  scale_y_continuous(breaks=seq(0,12,2), limits=c(0,12))
  
CairoPDF("diagnostic_yield_by_pheno.pdf")
p
dev.off()

# Check diagnostic yield per phenotype
s <- readLines("solved_platekey.txt")
q <- y %>% filter(grepl(paste(ka,collapse="|"), normalised_hpo_term, ignore.case=TRUE)) %>% select(plate_key)
r <- unique(gmc$plate_key[which(gmc$plate_key %in% q$plate_key)])
(length(unique(r))-length(unique(setdiff(r,s))))/length(unique(r))*100

both %>% filter(solved==1) %>% summarise(count=n_distinct(plate_key))
both %>% filter(solved==1 & fh==1) %>% summarise(count=n_distinct(plate_key))
both %>% filter(fh==1) %>% summarise(count=n_distinct(plate_key))

both %>% filter(solved==1 & age<18) %>% summarise(count=n_distinct(plate_key))

#Check previous genetic testing reports
query <- "SELECT rd.participant_id, gen.abnormal_cytogenetic_result, gen.abnormal_molecular_result, gen.method_of_test_id, gen.test_result_id, gen.test_scope
          FROM rare_diseases_invest_genetic_test_result AS gen
          LEFT JOIN rare_disease_analysis AS rd
          ON rd.participant_id = gen.participant_id
          WHERE rd.normalised_specific_disease = 'Congenital Anomaly of the Kidneys and Urinary Tract (CAKUT)'
          AND rd.participant_type = 'Proband'"

mysql <- labkey.executeSql(
  schemaName="lists",                                                 # Do not change this
  colNameOpt = "rname",                                               # Do not change this
  maxRows = 100000000,                                                # Do not change this
  folderPath="/main-programme/main-programme_v12_2021-05-06",          # This can be changed to different main programme releases
  sql = query                                                         # This can be changed to your query of choice
)

x <-data.frame(mysql)
test <- unique(x)


