#/bin/bash


CHROM=$1
VCF="results/vcfs/merged_calls.${CHROM}.vcf.gz"
REF="resources/ag1000g_WestAfrica_col_${CHROM}.vcf.gz"
MAP="resources/geneticMaps/${CHROM}.gmap"

while IFS="" read -r LINE || [ -n "$LINE" ];
do
printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
IRG=$(echo $LINE | cut -d" " -f3)
ORG=$(echo $LINE | cut -d" " -f4)
OUT="results/imputed_${CHROM}.${ID}.vcf"

workflow/scripts/GLIMPSE/static_bins/GLIMPSE_phase_static --input ${VCF} --reference ${REF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output ${OUT}

bcftools index ${OUT} 
done < "resources/chunks.${CHROM}.txt"