#!/usr/bin/env bash

#BSUB -q medium
#BSUB -P re_gecip_renal
#BSUB -o paintor_%J.stdout
#BSUB -e paintor_%J.stderr
#BSUB -cwd /re_gecip/renal/mchan/PUV/saige/saige_20210312/paintor
#BSUB -R "rusage[mem=10000]"

#Statistical fine mapping method that integrates functional genomic data with association strength from potentially multiple populations to prioritize variants for follow-up analysis.
#Input is a) summary association statistics (Zscores,with header +/- meta data - CHR POS RSID) , b) LD matrix (Pairwise Pearson correltations coefficients between each SNP), c) functional annotation matrix (with header)
# Outputs probability for a SNP to be causal to prioritize variants. Can model multiple causal variants at any risk locus. Leverage functional genomic data as a prior probability to improve prioritization.
# Also quantifies enrichment of causal variants within functional classes. Bayesian treatment of causal effect sizes.
# NB single space delimited hg19 only

module purge
module load bio/PAINTOR/v3.0
module load lang/Python/2.7.16-GCCcore-8.3.0

# Create locus file
chr=$1
start=$2
end=$3
trait=puv

grep "^chr${chr}" ../step2_v0.42.1_${trait}.SAIGE.gwas_concat_misHWEfiltered.txt | awk -v start="${start}" -v end="${end}" 'NR>1 {if ($2>start && $2<end) print $1, $2, $4, $5, $13, $10/$11}' > chr${chr}_locus_hg38.txt
sed -i '1i CHR pos A1 A2 P ZSCORE' chr${chr}_locus_hg38.txt

# Create LD matrix using 1000 genomes data as reference
for pop in EUR AFR SAS; do
python /resources/tools/manual_apps/software/bio/PAINTOR/v3.0-GCC-8.3.0/PAINTOR_V3.0/PAINTOR_Utilities/CalcLD_1KG_VCF.py \
--locus chr${chr}_locus_hg38_p0.05.txt \
--reference /public_data_resources/1000-genomes/20130502_GRCh38/ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz \
--map /public_data_resources/1000-genomes/20130502_GRCh37/integrated_call_samples_v3.20130502.ALL.panel \
--effect_allele A2 \
--alt_allele A1 \
--population ${pop} \
--Zhead ZSCORE \
--out_name chr${chr}_${pop}_p0.05 \
--position pos
mv chr${chr}_${pop}_p0.05.processed chr${chr}_${pop}_p0.05
done

# Create annotation file using Dnase, TFBS, CRE, gencode, CADD scores. Download from UCSC genome browser.
/public_data_resources/CADD/v1.6/GRCh38/whole_genome_SNVs_inclAnno.tsv.gz

for pop in EUR AFR SAS; do
python /resources/tools/manual_apps/software/bio/PAINTOR/v3.0-GCC-8.3.0/PAINTOR_V3.0/PAINTOR_Utilities/AnnotateLocus.py \
--input annotationpaths \
--locus chr${chr}_${pop}_p0.05 \
--out chr${chr}_${pop}_p0.05.annotations \
--chr CHR \
--pos pos
done

# input.files can be used to run on multiple loci at once
# Run with no annotations first and then run with each annotation individually, then pick top independent 4-5
PAINTOR -input input.files -Zhead ZSCORE -LDname ld -in ./ -out ./analysis -enumerate 2 -Gname Enrich.Base.chr12 -Lname BF.Base.chr12

while read ann; do PAINTOR -input input.files -Zhead ZSCORE -LDname ld -in ./ -out ./analysis -enumerate 2 -Gname Enrich.${ann} -Lname BF.${ann} -annotations ${ann}; done < annotations

PAINTOR -input input.files -Zhead ZSCORE -LDname ld -in ./ -out ./analysis -enumerate 3 -Gname Enrich.chr12_${pop} -Lname BF.chr12_${pop} -annotations dnase_hg38.bed,encode_cre.txt,tfclusters_hg38.bed,phastcons_hg38.bed,gencode_transcript_v29.bed,h1hesc_hic_hg38.bed
PAINTOR -input input.files -Zhead ZSCORE -LDname ld -in ./ -out ./analysis -enumerate 3 -Gname Enrich.chr6_${pop} -Lname BF.chr6_${pop} -annotations dnase_hg38.bed,encode_cre.txt,tfclusters_hg38.bed,phastcons_hg38.bed,gencode_transcript_v29.bed

