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

datadir=/gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/genomic_data_subset/PASS_MAF0.001/plinkfiles/plink_merged_by_chr
logs=${workdir}/logfiles
scripts=${workdir}/batchscripts

mkdir -p $logs
mkdir -p $scripts

files=$(ls ${datadir}/*.bim)

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
module load bio/magma/1.09a

magma --annotate --snp-loc $SNPLOC_FILE --gene-loc $GENELOC --out $ANNOT_PREFIX
magma --bfile [DATA] --gene-annot [ANNOT_PREFIX].genes.annot --out [GENE_PREFIX]


    bsub < ${scripts}/script_${trait}_${chunk}.sh

done
SNPLOC_FILE=
GENELOC=/re_gecip/renal/mchan/CAKUT/lists/NCBI38.gene.loc
ANNOT_PREFIX=cakut

magma --annotate --snp-loc $SNPLOC_FILE --gene-loc $GENELOC --out $ANNOT_PREFIX

magma --bfile [DATA] --gene-annot [ANNOT_PREFIX].genes.annot --out [GENE_PREFIX]

magma --bfile [REFDATA] --pval [PVAL_FILE] N=[N] --gene-annot [ANNOT_PREFIX].genes.annot --out [GENE_PREFIX]

magma --gene-results [GENE_PREFIX].genes.raw --set-annot [SET_FILE] --out [GS_PREFIX]

magma --gene-results [GENE_PREFIX].genes.raw --set-annot [SET_FILE] --out [GS_PREFIX]
