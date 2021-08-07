#/bin/bash

CHROM=$1
VCF=$2
REF=$3
MAP=$4
CHUNK=$5
THREADS=$6
DATASET=$7

while IFS="" read -r LINE || [ -n "$LINE" ];
do
printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
IRG=$(echo $LINE | cut -d" " -f3)
ORG=$(echo $LINE | cut -d" " -f4)
OUT="results/${dataset}VCFs/imputed.${CHROM}.${ID}.vcf"

workflow/scripts/GLIMPSE/static_bins/GLIMPSE_phase_static --input ${VCF} --reference ${REF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output ${OUT} --thread ${THREADS}

bcftools index ${OUT}
done < $CHUNK

find results/${dataset}VCFs/ -name "imputed.${CHROM}*.vcf" | sort > results/${dataset}VCFs/${CHROM}.vcf.list

echo "Imputation complete..."

LST=results/${dataset}VCFs/${CHROM}.vcf.list
OUT=results/${dataset}VCFs/imputed.${CHROM}.vcf.gz

workflow/scripts/GLIMPSE/static_bins/GLIMPSE_ligate_static --input ${LST} --output ${OUT} --thread ${THREADS}

echo "Ligation complete..."
