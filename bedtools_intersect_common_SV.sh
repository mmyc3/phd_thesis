#!/usr/bin/env bash

#BSUB -q medium
#BSUB -P re_gecip_renal
#BSUB -o sv_common_%J_%I.stdout
#BSUB -e sv_common_%J_%I.stderr
#BSUB -J "sv[1-4]"
#BSUB -cwd /re_gecip/renal/mchan/CAKUT/SV/data
#BSUB -R "rusage[mem=10000] span[hosts=1]"

module load bio/BEDTools/2.27.1-foss-2018b

# This script loops over individual SV VCFs to remove SVs that overlap with a set of common SVs from dbVar (> 1%) and an in-house control set
# Split separately into SV type as called by CANVAS and MANTA - CNV, INV, DEL, DUP

SV=$(sed -n ${LSB_JOBINDEX}p sv_list.txt)
COHORT=cakut_ancestry_matched_controls
GENE=gencode
SCRATCH="/re_scratch/re_gecip/renal/mchan"

VCFS="/re_gecip/renal/mchan/CAKUT/SV/data/${COHORT}_${GENE}_vcf.txt"
COMMON="/re_gecip/renal/mchan/CAKUT/SV/data/dbvar_${SV}_hg38_final.bed"
CANCER="/re_gecip/renal/mchan/CAKUT/SV/cancer/cancer_${SV}.bed"
OUTPATH="/re_gecip/renal/mchan/CAKUT/SV/cakut/analysis"

while read -r vcf; do
	vfile="${vcf##*/}"
        OUTFILE="${OUTPATH}/${vfile}.${SV}_new.vcf"
	grep "^#" /re_gecip/renal/mchan/CAKUT/SV/data/${vcf} > ${SCRATCH}/${vfile}.${GENE}.header
	grep ${SV} /re_gecip/renal/mchan/CAKUT/SV/data/${vcf} > ${SCRATCH}/${SV}.${vfile}
	cat ${SCRATCH}/${vfile}.${GENE}.header ${SCRATCH}/${SV}.${vfile} > ${SCRATCH}/${vfile}.${SV}.vcf
	bedtools intersect -a ${SCRATCH}/${vfile}.${SV}.vcf \
	-b ${CANCER} \
	-v \
	-r \
	-f 0.7 \
	-header > ${OUTFILE}
done < ${VCFS}

cd $OUTPATH
ls *${GENE}.vcf.${SV}_new.vcf > ${GENE}_${SV}.txt
