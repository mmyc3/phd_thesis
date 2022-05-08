##### This script uses standalone LocusZoom to plot GWAS association results #####
# NB there are no recombination maps for GRCh38 #####
##refFlat is deafult gene table, can change to gencode
## SNP names should be rsIDs or CHR:POS format (i.e. need to remove ref:alt from ID)

#!/usr/bin/env bash

#BSUB -q short
#BSUB -P re_gecip_renal
#BSUB -o zoom_%J.stdout
#BSUB -e zoom_%J.stderr
#BSUB -cwd /re_gecip/renal/mchan/CAKUT/saige/bladderexstrophy
#BSUB -R "rusage[mem=10000]"

module purge
module load lang/R/3.6.2-foss-2019b
module load bio/HTSlib/1.9-foss-2018b
module load bio/PLINK/1.9b_4.1-x86_64
module load bio/locuszoom/1.4-R_3.6.2-HTSlib_1.9

COHORT=cakut
IN=${COHORT}_gwas_newIDs
SNP=$(sed -n ${LSB_JOBINDEX}p zoom_snps.txt)

locuszoom --metal ${IN} --build hg38 --pop EUR --source 1000G_Nov2014 --refsnp ${SNP} --flank 500kb --prefix ${COHORT}

#--bed-tracks ${BED}
#--chr ${CHR} --start ${POS_FROM} --end ${POS_TO}
#--refgene ${GENE} --flank 250kb \
