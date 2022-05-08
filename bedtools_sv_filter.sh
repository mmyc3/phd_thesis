#############################################################################################################################

#USE Bedtools to filter individual SV vcf files #################################################################

# This script subsets VCFs by coordinate and PASS filter 

# INPUT:
        # 1: Tab-delimited list of chromosomal coordinates (regions.txt)
        # 2. A list of full-paths to VCF files of interest (vcf_list.txt)

#! usr/bin/env bash

#BSUB -q medium
#BSUB -P re_gecip_renal
#BSUB -o sv_%J.stdout
#BSUB -e sv_%J.stderr
#BSUB -cwd /re_gecip/renal/mchan/PUV/data/
#BSUB -R "rusage[mem=10000]"

module load bio/BEDTools/2.27.1-foss-2018b

VCFS="/re_gecip/renal/mchan/PUV/data/${COHORT}_vcf_list.txt"
REGIONS="/re_gecip/renal/mchan/PUV/data/${GENE}.bed"
OUTPATH="/re_gecip/renal/mchan/PUV/analysis/"

while read -r vcf; do
	vfile="${vcf##*/}"
        OUTFILE="${OUTPATH}/${vfile}.${COHORT}_${GENE}.vcf"
	bedtools intersect -a ${vcf} \
	-b ${REGIONS} \
	-u \
	-header > ${OUTFILE}
done < ${VCFS}

cd $OUTPATH
ls *.${COHORT}_${GENE}.vcf > ${COHORT}_vcf.txt

while read -r vcf; do
	vfile="${vcf##*/}"
        OUTFILE="${OUTPATH}/${vfile}.filtered.vcf"
	bcftools view ${vcf} \
	-f PASS,. \
	-o ${OUTFILE}
done < ${COHORT}_vcf.txt

rm *${COHORT}_${GENE}.vcf
ls *.${COHORT}_${GENE}.vcf.filtered.vcf > ${COHORT}_vcf.txt
sed -i 's/^/\/re_gecip\/renal\/mchan\/PUV\/analysis\//g' ${COHORT}_vcf.txt
