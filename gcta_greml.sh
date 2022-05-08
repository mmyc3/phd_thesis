#!/usr/bin/env bash

#BSUB -q short
#BSUB -P re_gecip_renal
#BSUB -o gcta_%J.stdout
#BSUB -e gcta_%J.stderr
#BSUB -cwd /re_gecip/renal/mchan/CAKUT/saige/eur/h2
#BSUB -R rusage[mem=10000] span[hosts=1] -M 10000 -n2

# GREML-LDMS estimates heritability using WGS data (Yang Nat Genet 2015)
# Corrects for LD bias in estimated SNP-based heritability
# Unbiased regardless of MAF or LD of underlying causal variants

# 1. Calculate segment based LD score

module load bio/GCTA/1.93.1_beta
module load bio/PLINK/2.00-devel-20200409-x86_64

workdir=/re_gecip/renal/mchan/CAKUT/saige/eur/h2
datadir=/gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/genomic_data_subset/PASS_MAF0.001/plinkfiles
phenofile="/re_gecip/renal/mchan/CAKUT/pheno/cakut_eur_ancestry_matched_unrelated_final.ped"
trait=cakut
scratch="/re_scratch/re_gecip/renal/mchan"
logs=${workdir}/logfiles
scripts=${workdir}/batchscripts

mkdir -p $logs
mkdir -p $scripts

files=$(ls ${datadir}/*.bim)

cd $workdir

for file in ${files}; do

    chunk=$(basename $file | sed 's/gel_mainProgramme_aggV2_//g' | sed 's/_filteredPASS.bim//g')
    CHR=$(echo $chunk | cut -f 1 -d "_")

    echo -e '#!/usr/bin/env bash' >${scripts}/script_${trait}_${chunk}.sh

    printf "#BSUB -q medium
#BSUB -P re_gecip_renal
#BSUB -J ${trait}_${chunk}_plink
#BSUB -o $logs/${trait}_${chunk}_plink.out
#BSUB -e $logs/${trait}_${chunk}_plink.err
#BSUB -R rusage[mem=10000] span[hosts=1] -M 10000 -n2

module purge
module load bio/PLINK/1.9b_4.1-x86_64
module load bio/GCTA/1.93.1_beta

plink --bfile $datadir/gel_mainProgramme_aggV2_${chunk}_filteredPASS --keep samples.txt --hwe 0.000001 --geno 0.01 --allow-no-sex --make-bed --memory 10000 --out $scratch/${chunk}
gcta64 --bfile $scratch/${chunk} --ld-score-region 200 --out $scratch/${chunk}" >> ${scripts}/script_${trait}_${chunk}.sh

    bsub < ${scripts}/script_${trait}_${chunk}.sh

done

for i in {1..22}; do find $scratch -name "chr${i}_*.score.ld" | xargs -n 1 tail -n +2 | sort -k1,1g > $scratch/chr${i}_all; cat header $scratch/chr${i}_all > chr${i}.score.ld; done
cat chr{[0-9],[0-9][0-9]}.score.ld >> all.score.ld

# # 2. Stratify SNPs based on LD score in R
Rscript --vanilla /re_gecip/renal/mchan/PUV/scripts/gcta_greml.R

for i in {1..22}; do plink2 --bfile /gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/genomic_data_subset/PASS_MAF0.001/plinkfiles/plink_merged_by_chr/gel_mainProgramme_aggV2_chr${i}_filteredPASS --keep plink_samples.txt --geno 0.01 --hwe 0.000001 --make-bed --out $scratch/chr${i}; done
ls $scratch/chr{[0-9],[0-9][0-9]}.bim | sed 's/.bim//g' | sort -V > plink_files.txt

# # 3. Compute GRMs using stratified SNPs based on LD quartile and MAF
bsub < /re_gecip/renal/mchan/scripts/gcta_grm.sh

# # 4. Perform REML analysis using multiple GRMs
awk 'BEGIN{OFS="\t"}{print $1,$2,$6}' $phenofile > $trait.pheno
awk 'BEGIN{OFS="\t"}{print $1,$2,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' $phenofile > $trait.pcs

ls *.grm.bin | sed 's/.grm.bin//g' > multi_grm.txt
gcta64 --reml --mgrm multi_grm.txt --pheno $trait.pheno --qcovar $trait.pcs --out puv --prevalence 0.002 --thread-num 10
