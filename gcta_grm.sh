#!/usr/bin/env bash

#BSUB -q medium
#BSUB -P re_gecip_renal
#BSUB -o gcta_%J_${i}.stdout
#BSUB -e gcta_%J_${i}.stderr
#BSUB -J "gcta[1-4]"
#BSUB -cwd /re_gecip/renal/mchan/CAKUT/h2/eur
#BSUB -R rusage[mem=10240] span[hosts=1] -M 10240 -n12

module load bio/GCTA/1.93.1_beta

scratch="/re_scratch/re_gecip/renal/mchan"
i=$(sed -n ${LSB_JOBINDEX}p group.txt)

gcta64 --mbfile plink_files.txt --extract snp_group$i.txt --make-grm-part 3 1 --maf 0.4 --max-maf 0.5 --out $scratch/group${i}_maf0.4 --thread-num 5
gcta64 --mbfile plink_files.txt --extract snp_group$i.txt --make-grm-part 3 2 --maf 0.4 --max-maf 0.5 --out $scratch/group${i}_maf0.4 --thread-num 5
gcta64 --mbfile plink_files.txt --extract snp_group$i.txt --make-grm-part 3 3 --maf 0.4 --max-maf 0.5 --out $scratch/group${i}_maf0.4 --thread-num 5

gcta64 --mbfile plink_files.txt --extract snp_group$i.txt --make-grm-part 3 1 --maf 0.3 --max-maf 0.4 --out $scratch/group${i}_maf0.3 --thread-num 5
gcta64 --mbfile plink_files.txt --extract snp_group$i.txt --make-grm-part 3 2 --maf 0.3 --max-maf 0.4 --out $scratch/group${i}_maf0.3 --thread-num 5
gcta64 --mbfile plink_files.txt --extract snp_group$i.txt --make-grm-part 3 3 --maf 0.3 --max-maf 0.4 --out $scratch/group${i}_maf0.3 --thread-num 5

gcta64 --mbfile plink_files.txt --extract snp_group$i.txt --make-grm-part 3 1 --maf 0.2 --max-maf 0.3 --out $scratch/group${i}_maf0.2 --thread-num 5
gcta64 --mbfile plink_files.txt --extract snp_group$i.txt --make-grm-part 3 2 --maf 0.2 --max-maf 0.3 --out $scratch/group${i}_maf0.2 --thread-num 5
gcta64 --mbfile plink_files.txt --extract snp_group$i.txt --make-grm-part 3 3 --maf 0.2 --max-maf 0.3 --out $scratch/group${i}_maf0.2 --thread-num 5

gcta64 --mbfile plink_files.txt --extract snp_group$i.txt --make-grm-part 3 1 --maf 0.1 --max-maf 0.2 --out $scratch/group${i}_maf0.1 --thread-num 5
gcta64 --mbfile plink_files.txt --extract snp_group$i.txt --make-grm-part 3 2 --maf 0.1 --max-maf 0.2 --out $scratch/group${i}_maf0.1 --thread-num 5
gcta64 --mbfile plink_files.txt --extract snp_group$i.txt --make-grm-part 3 3 --maf 0.1 --max-maf 0.2 --out $scratch/group${i}_maf0.1 --thread-num 5

gcta64 --mbfile plink_files.txt --extract snp_group$i.txt --make-grm-part 3 1 --maf 0.05 --max-maf 0.1 --out $scratch/group${i}_maf0.05 --thread-num 5
gcta64 --mbfile plink_files.txt --extract snp_group$i.txt --make-grm-part 3 2 --maf 0.05 --max-maf 0.1 --out $scratch/group${i}_maf0.05 --thread-num 5
gcta64 --mbfile plink_files.txt --extract snp_group$i.txt --make-grm-part 3 3 --maf 0.05 --max-maf 0.1 --out $scratch/group${i}_maf0.05 --thread-num 5

gcta64 --mbfile plink_files.txt --extract snp_group$i.txt --make-grm-part 3 1 --maf 0.01 --max-maf 0.05 --out $scratch/group${i}_maf0.01 --thread-num 5
gcta64 --mbfile plink_files.txt --extract snp_group$i.txt --make-grm-part 3 2 --maf 0.01 --max-maf 0.05 --out $scratch/group${i}_maf0.01 --thread-num 5
gcta64 --mbfile plink_files.txt --extract snp_group$i.txt --make-grm-part 3 3 --maf 0.01 --max-maf 0.05 --out $scratch/group${i}_maf0.01 --thread-num 5
#
gcta64 --mbfile plink_files.txt --extract snp_group$i.txt --make-grm-part 3 1 --maf 0.001 --max-maf 0.01 --out $scratch/group${i}_maf0.001 --thread-num 5
gcta64 --mbfile plink_files.txt --extract snp_group$i.txt --make-grm-part 3 2 --maf 0.001 --max-maf 0.01 --out $scratch/group${i}_maf0.001 --thread-num 5
gcta64 --mbfile plink_files.txt --extract snp_group$i.txt --make-grm-part 3 3 --maf 0.001 --max-maf 0.01 --out $scratch/group${i}_maf0.001 --thread-num 5

for maf in 0.001 0.01 0.05 0.1 0.2 0.3 0.4; do
cat ${scratch}/group${i}_maf${maf}.part_3_*.grm.id > group${i}_maf${maf}.grm.id
cat ${scratch}/group${i}_maf${maf}.part_3_*.grm.bin > group${i}_maf${maf}.grm.bin
cat ${scratch}/group${i}_maf${maf}.part_3_*.grm.N.bin > group${i}_maf${maf}.grm.N.bin
done
