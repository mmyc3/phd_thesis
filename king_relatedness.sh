#!/usr/bin/env bash

#BSUB -q short
#BSUB -P re_gecip_renal
#BSUB -o king_%J_%I.stdout
#BSUB -e king_%J_%I.stderr
#BSUB -J "king[1-5]"
#BSUB -cwd /re_gecip/renal/mchan/CAKUT/pheno
#BSUB -R rusage[mem=10000]

module load bio/PLINK/2.00-devel-20200409-x86_64
module load bio/KING/2.2.4
module load lang/Python/3.7.4-GCCcore-8.3.0

# This script uses the KING-robust algorithm in PLINK 2 to create a set of unrelated (greater than second degree relationships) cases and controls. 
# The input files are a set of ~120,000 high-quality LD-pruned variants with MAF > 1%
# The kinship coefficients are recalculated on the combined case-control cohort and a custom python script used to preferentially keep cases over controls. 
# PCA is then run on the final unrelated case-control set to generate 10 principal components.

CASES=cakut
$(sed -n ${LSB_JOBINDEX}p pheno_list.txt)

# unrelated controls
plink2 \
--bfile /gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/additional_data/HQ_SNPs/MAF1/GELautosomes_LD_pruned_1kgp3Intersect_maf0.01_mpv10 \
--keep controls.txt \
--king-cutoff 0.0884 \
--out controls

# unrelated cases
plink2 \
--bfile /gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/additional_data/HQ_SNPs/MAF1/GELautosomes_LD_pruned_1kgp3Intersect_maf0.01_mpv10 \
--keep $CASES.txt \
--king-cutoff 0.0884 \
--out $CASES

# combine cases and controls 
cat controls.king.cutoff.in.id $CASES.king.cutoff.in.id | grep -v "FID" > ${CASES}_control.txt

# extract case-control samples 
plink2 --bfile /gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/additional_data/HQ_SNPs/MAF1/GELautosomes_LD_pruned_1kgp3Intersect_maf0.01_mpv10 \
--keep ${CASES}_control.txt \
--make-bed \
--out ${CASES}_control

# recalculate kinship coefficients using KING (Manichaikul et al, 2010)
king -b ${CASES}_control.bed --related --degree 2 --prefix ${CASES}

# if cases and controls are related, preferentially keep cases
/re_gecip/renal/mchan/scripts/remove-related-controls.py ${CASES}.king.cutoff.in.id controls.king.cutoff.in.id king.kin0 keep_${CASES}.txt
grep -v FID keep_${CASES}.txt > ${CASES}_control_unrelated.txt

# run PCA on the final unrelated case-control cohort 
plink2 --bfile /gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/additional_data/HQ_SNPs/MAF1/GELautosomes_LD_pruned_1kgp3Intersect_maf0.01_mpv10 \
--keep case_control_unrelated.txt \
--pca 10
