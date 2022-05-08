#!/usr/bin/env bash

#BSUB -q short
#BSUB -P re_gecip_renal
#BSUB -o magma_%J_%I.stdout
#BSUB -e magma_%J_%I.stderr
#BSUB -J "magma[1-23]"
#BSUB -cwd /re_gecip/renal/mchan/CAKUT/saige/magma
#BSUB -R rusage[mem=10000] span[hosts=1] -M 10000 -n2

trait=cakut
scratch=/re_scratch/re_gecip/renal/mchan/
phenofile=/re_gecip/renal/mchan/CAKUT/pheno/cakut_ancestry_matched_controls_unrelated_aggv2_pcs.ped
datadir=/gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/genomic_data_subset/PASS_MAF0.001/plinkfiles/plink_merged_by_chr
pheno=pheno
covariates=covar

i=$(sed -n ${LSB_JOBINDEX}p chr_list.txt)

module purge
module load bio/magma/1.09a

magma --bfile $datadir/gel_mainProgramme_aggV2_chr${i}_filteredPASS --gene-annot $trait.genes.annot --gene-model multi=all --covar file=$covariates --pheno file=$pheno use=PHENO --gene-settings adap-permp --burden 0.01 --out $trait.$i
