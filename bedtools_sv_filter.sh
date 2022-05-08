#############################################################################################################################

#USE Bedtools to filter individual SV vcf files #################################################################

# This script subsets VCFs by coordinate and PASS filter 

# INPUT:
        # 1: Tab-delimited list of chromosomal coordinates (regions.txt)
        # 2. A list of full-paths to VCF files of interest (vcf_list.txt)
	
#! /usr/bin/env bash

#BSUB -q long
#BSUB -P re_gecip_renal
#BSUB -o sv_%J.stdout
#BSUB -e sv_%J.stderr
#BSUB -cwd /re_gecip/renal/mchan/CAKUT/SV/cakut
#BSUB -R "rusage[mem=10000] span[hosts=1]"

module load bio/BEDTools/2.27.1-foss-2018b
module load bio/BCFtools/1.11-GCC-8.3.0

COHORT=cakut_ancestry_matched_controls
GENE=gencode

VCFS="/re_gecip/renal/mchan/CAKUT/SV/data/cakut_control_vcf_list_v10.txt"
EXONS="/re_gecip/renal/mchan/CAKUT/SV/data/gencode_exons_v29.bed"
CRE="/re_gecip/renal/mchan/CAKUT/lists/encode_cre.txt"
OUTPATH="/re_gecip/renal/mchan/CAKUT/SV/cakut/data"

while read -r vcf; do
	vfile="${vcf##*/}"
        OUTFILE="${OUTPATH}/${vfile}.${GENE}.vcf"
	bedtools intersect -a ${vcf} -b ${EXONS} -u -header | \
	bcftools view -f PASS -o ${OUTFILE}
done < ${VCFS}

cd ${OUTPATH}
ls *.${GENE}.vcf > ${COHORT}_${GENE}_vcf.txt


