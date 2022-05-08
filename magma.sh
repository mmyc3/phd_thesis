#!/usr/bin/env bash

#BSUB -q short
#BSUB -P re_gecip_renal
#BSUB -o magma_%J.stdout
#BSUB -e magma_%J.stderr
#BSUB -cwd /re_gecip/renal/mchan/CAKUT/saige/magma
#BSUB -R rusage[mem=10000] span[hosts=1] -M 10000 -n2

# MAGMA carries out gene-level analysis using GWAS data or raw SNP data
#first, an annotation step to map SNPs onto genes;
#second, a gene analysis step to compute gene p-values;
#third, a gene-level analysis step: either a generalized gene-set analysis, a gene property analysis, or both.

trait=cakut
scratch=/re_scratch/re_gecip/renal/mchan/
phenofile=/re_gecip/renal/mchan/CAKUT/pheno/cakut_ancestry_matched_controls_unrelated_aggv2_pcs.ped
datadir=/gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/genomic_data_subset/PASS_MAF0.001/plinkfiles/plink_merged_by_chr
GENELOC=/re_gecip/renal/mchan/CAKUT/lists/NCBI38.gene.loc
hallmark=/re_gecip/renal/mchan/CAKUT/lists/hallmark_v7.4.txt
go=/re_gecip/renal/mchan/CAKUT/lists/c5.go.v7.4.symbols.gmt.txt
reg=/re_gecip/renal/mchan/CAKUT/lists/c3.all.v7.4.symbols.gmt.txt
hpo=/re_gecip/renal/mchan/CAKUT/lists/c5.hpo.v7.4.symbols.gmt.txt
cp=/re_gecip/renal/mchan/CAKUT/lists/c2.cp.v7.4.symbols.gmt.txt

module purge
module load bio/magma/1.09a

# Annotate files to genes (NB doesn't recognise chr prefix)

awk 'BEGIN{OFS="\t"} {print $6, $2, $3, $4}' $GENELOC > gene_loc.txt

for i in {1..22} X; do
  sed 's/^chr//g' $datadir/gel_mainProgramme_aggV2_chr${i}_filteredPASS.bim > $scratch/chr${i}.bim
  magma --annotate --snp-loc $scratch/chr${i}.bim --gene-loc gene_loc.txt --out ${trait}_chr${i}
  rm $scratch/chr${i}.bim
done

# Merge all annotation files into one

head -2 ${trait}_chr22.genes.annot > $trait.genes.annot

for i in {1..22} X; do
  tail -n+3 ${trait}_chr${i}.genes.annot >> $trait.genes.annot
done


# Gene analysis
# Uses 1/2 for control/case
# Default burden test for AF < 1% - weighted by inverse AF. Combined with common variants. Can add own weights.
# max allowed SNP missingness is 0.05
# differential missingness removed is p < 1e-6
# can calculate empirical p values via permutation - uses adaptive permutation, varying no of permutations by gene
# gene-model default is linear regression using PCs - regresses the phenotype on principal components derived from SNPs in gene. Can combine all into aggregate p value

awk '{print $1,$2,$6}' $phenofile > pheno
sed -i 's/1$/2/g' pheno
sed -i 's/0$/1/g' pheno
pheno=pheno

awk '{print $1,$2,$5,$7, $8, $9, $10, $11, $12, $13, $14, $15, $16}' $phenofile > covar
covariates=covar

bsub < /re_gecip/renal/mchan/scripts/magma_gene_analysis_per_chr.sh

for i in {1..22} X; do
  awk '{if($15<0.05) print $0}' ${trait}.${i}.genes.out >> $trait.genes.sig
done

# Gene-set analysis
head -2 ${trait}.22.genes.raw > $trait.genes.raw

for i in {1..22} X; do
  tail -n+3 ${trait}.${i}.genes.raw >> $trait.genes.raw
done

magma --gene-results $trait.genes.raw --set-annot $hallmark --out $trait.hallmark
magma --gene-results $trait.genes.raw --set-annot $go --out $trait.go
magma --gene-results $trait.genes.raw --set-annot $cp --out $trait.cp
magma --gene-results $trait.genes.raw --set-annot $hpo --out $trait.hpo
magma --gene-results $trait.genes.raw --set-annot $reg --out $trait.reg
