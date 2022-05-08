#!/usr/bin/env bash

#BSUB -q short
#BSUB -P re_gecip_renal
#BSUB -o epistasis_%J.stdout
#BSUB -e epistasis_%J.stderr
#BSUB -cwd /re_gecip/renal/mchan/PUV/saige/epistasis
#BSUB -R "rusage[mem=10000] span[hosts=1]"

module load bio/PLINK/1.9b_4.1-x86_64

SET="variant_set.txt"
PHENO="/re_gecip/renal/mchan/PUV/pca/puv_ancestry_matched_controls_unrelated_aggv2_final.ped"

plink --bfile ../fine_map/chr6_12 --pheno ${PHENO} --pheno-name PHENO --1 --allow-no-sex --epistasis set-by-set --epi1 1 --out puv --set ${SET}
