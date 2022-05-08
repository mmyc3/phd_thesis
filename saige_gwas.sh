#!/usr/bin/env bash

#BSUB -q medium
#BSUB -P re_gecip_renal
#BSUB -o gwas_%J.stdout
#BSUB -e gwas_%J.stderr
#BSUB -cwd /re_gecip/renal/mchan/PUV/saige/
#BSUB -R "rusage[mem=10000]"

module purge
module load lang/R/4.0.2-foss-2019b
module load singularity/3.2.1

trait=puv

plink_folder_for_step1=/gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/additional_data/HQ_SNPs/
plink_files_for_step1=${plink_folder_for_step1}/GELautosomes_LD_pruned_1kgp3Intersect_maf0.05_mpv10
outputdir=SAIGE_step1
workdir=/re_gecip/renal/mchan/PUV/saige/
datadir=/gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/genomic_data_subset/PASS_MAF0.001
step1_dir=SAIGE_step1
males_sample_file=${workdir}/pheno_maleslist.txt
phenofile="/re_gecip/renal/mchan/PUV/pca/puv_ancestry_matched_unrelated_aggv2_pcs_final.ped"
phenotypeCol="PHENO"
threads=1
SAIGE_version="0.42.1"

mkdir -p $workdir
cd $workdir
mkdir -p $outputdir

singularity exec -B /gel_data_resources /gel_data_resources/containers/singularity_containers/saige_0.42.1.sif step1_fitNULLGLMM.R \
    --plinkFile=${plink_files_for_step1} \
    --phenoFile=${phenofile} \
    --phenoCol=${phenotypeCol} \
    --sampleIDColinphenoFile=IID \
    --traitType=binary \
    --outputPrefix=${outputdir}/step1_v${SAIGE_version}_${trait} \
    --outputPrefix_varRatio=${outputdir}/step1_v${SAIGE_version}_${trait} \
    --nThreads=${threads} \
    --LOCO=FALSE \
    --IsOverwriteVarianceRatioFile=TRUE \
    --covarColList=SEX,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10

awk '{if($5==1) print $2}' $phenofile > pheno_maleslist.txt

outputdir=$workdir/SAIGE_step2
files=$(ls ${datadir}/*.vcf.gz)
cd $workdir

### jobs folders
logs=${workdir}/logfiles/SAIGE_step2
scripts=${workdir}/batchscripts/SAIGE_step2

mkdir -p $logs
mkdir -p $scripts
mkdir -p $outputdir

for file in ${files}; do

    chunk=$(basename $file | sed 's/gel_mainProgramme_aggV2_//g' | sed 's/_filteredPASS.vcf.gz//g')
    CHR=$(echo $chunk | cut -f 1 -d "_")

    echo -e '#!/usr/bin/env bash' >${scripts}/script_${trait}_${chunk}.sh

    printf "#BSUB -q short
#BSUB -P re_gecip_renal
#BSUB -J ${trait}_${chunk}_SAIGE_step2
#BSUB -o $logs/${trait}_${chunk}_SAIGE_step2.out
#BSUB -e $logs/${trait}_${chunk}_SAIGE_step2.err

module purge
module load singularity/3.2.1

singularity exec -B ${datadir} /gel_data_resources/containers/singularity_containers/saige_0.41.sif step2_SPAtests.R \
--vcfFile=${file} \
--vcfFileIndex=${file}.tbi \
--vcfField=GT \
--chrom=${CHR} \
--minMAC=20 \
--GMMATmodelFile=${step1_dir}/step1_v${SAIGE_version}_${trait}.rda \
--varianceRatioFile=${step1_dir}/step1_v${SAIGE_version}_${trait}.varianceRatio.txt \
--SAIGEOutputFile=${outputdir}/step2_v${SAIGE_version}_${trait}.${chunk}.SAIGE.gwas.txt \
--numLinesOutput=2 \
--IsOutputAFinCaseCtrl=TRUE \
--IsDropMissingDosages=TRUE \
--IsOutputHetHomCountsinCaseCtrl=TRUE \
--IsOutputNinCaseCtrl=TRUE" >> ${scripts}/script_${trait}_${chunk}.sh

if [ ${CHR} = "chrX" ]; then
printf " --sampleFile_male=${males_sample_file} \
--X_PARregion=10001-2781479,155701383-156030895 \
--is_rewrite_XnonPAR_forMales=TRUE
" >> ${scripts}/script_${trait}_${chunk}.sh
fi
    bsub -R 'rusage[mem=10000]' < ${scripts}/script_${trait}_${chunk}.sh

done

### STEP 2 use PLINK for QC filtering 

inputdir=$workdir/SAIGE_step2
outputdir=$workdir

cd $workdir
files=$(ls ${inputdir}/step2_v${SAIGE_version}_${trait}.*.SAIGE.gwas.txt)
file1=$(echo $files|cut -f 1 -d ' ')
head -n 1 $file1 > ${outputdir}/step2_v${SAIGE_version}_${trait}.SAIGE.gwas_concat.txt

for file in ${files}; do

cat $file|grep -v SNPID >> ${outputdir}/step2_v${SAIGE_version}_${trait}.SAIGE.gwas_concat.txt

done

datadir=/gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/genomic_data_subset/PASS_MAF0.001/plinkfiles
outputdir=$workdir/siteQC

cd $workdir

thresm=0.00001
threshwe=0.000001

files=$(ls ${datadir}/*.bed)

mkdir -p $outputdir
mkdir -p $outputdir/logfile
mkdir -p $outputdir/batchscripts

for file in ${files}; do

    file0=$(basename $file | sed 's/.bed//g')
    outfile=$outputdir/$file0

    chunk=$(basename $file | sed 's/gel_mainProgramme_aggV2_//g' | sed 's/_filteredPASS.bed//g')
    CHR=$(echo $chunk | cut -f 1 -d "_")

    echo -e '#!/usr/bin/env bash' >${outputdir}/batchscripts/${file0}_siteQC.sh
    echo "#BSUB -q short
#BSUB -P re_gecip_renal
#BSUB -o $outputdir/logfile/${file0}_siteQC.out
#BSUB -e $outputdir/logfile/${file0}_siteQC.err

module load bio/PLINK/1.9b_4.1-x86_64

plink --bfile ${datadir}/${file0} --out ${outfile} --1 --allow-no-sex --keep-allele-order --geno 0.01 --test-missing --pheno ${phenofile} --pheno-name ${phenotypeCol} --output-chr chrM --memory 5000

awk -v thresm=$thresm '\$5 < thresm {print}' ${outfile}.missing > ${outfile}.missing_FAIL
" >>${outputdir}/batchscripts/${file0}_siteQC.sh

    if [ ${CHR} = "chrX" ]; then
        echo "
plink --bfile ${datadir}/${file0} --out ${outfile}.misHWE --1 --allow-no-sex --keep-allele-order --geno 0.01 --hwe ${threshwe} --exclude ${outfile}.missing_FAIL --pheno ${phenofile} --pheno-name ${phenotypeCol} --make-just-bim --memory 5000 --update-sex ${phenofile} 3 --filter-females --output-chr chrM

" >>${outputdir}/batchscripts/${file0}_siteQC.sh

    else
        echo "
plink --bfile ${datadir}/${file0} --out ${outfile}.misHWE --1 --allow-no-sex --keep-allele-order --geno 0.01 --hwe ${threshwe} --exclude ${outfile}.missing_FAIL --pheno ${phenofile} --pheno-name ${phenotypeCol} --make-just-bim --output-chr chrM --memory 5000

" >>${outputdir}/batchscripts/${file0}_siteQC.sh

    fi

    bsub -R 'rusage[mem=10000]' <${outputdir}/batchscripts/${file0}_siteQC.sh

done

### STEP 3 - concat filtered results and plot Manhattan and Q-Q 

inputdir=$workdir/siteQC
outputdir=$workdir
outputtag=aggV2_${trait}_concat

cd $workdir
files=$(ls ${inputdir}/*.misHWE.bim)
>${outputdir}/${outputtag}.misHWE_PASS
for file in ${files}; do

cat $file >> ${outputdir}/${outputtag}.misHWE_PASS

done

head -n 1 step2_v0.42.1_${trait}.SAIGE.gwas_concat.txt > step2_v0.42.1_${trait}.SAIGE.gwas_concat_misHWEfiltered.txt
awk 'NR==FNR{c[$1$4$6$5]++;next}; c[$1$2$4$5] > 0' aggV2_${trait}_concat.misHWE_PASS step2_v0.42.1_${trait}.SAIGE.gwas_concat.txt >> step2_v0.42.1_${trait}.SAIGE.gwas_concat_misHWEfiltered.txt

module load lang/R/3.6.2-foss-2019b

Rscript --vanilla /re_gecip/renal/mchan/PUV/scripts/saige_manhattan_qq.R $trait

