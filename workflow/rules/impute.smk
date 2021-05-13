# An example collection of Snakemake rules imported in the main Snakefile.

rule glimpseChunk:
    input:
        sites = "resources/ag1000g_WestAfrica_col_{chrom}.sites.vcf.gz",
    output:
        chunks = "resources/chunks.{chrom}.txt"
    log:
        "logs/glimpseChunk/{chrom}.log"
    shell:
        """
        workflow/scripts/GLIMPSE/static_bins/GLIMPSE_chunk_static --input {input.sites} --region {wildcards.chrom} --window-size 2000000 \
        --buffer-size 200000 --output {output} 2> {log}
        """


rule glimpseImpute:
    input:
        vcf = "results/vcfs/merged_calls.{chrom}.vcf.gz",
        HapPanel = "resources/ag1000g_WestAfrica_col_{chrom}.vcf.gz",
        genMap = "resources/geneticMaps/{chrom}.gmap",
        chunks = "resources/chunks.{chrom}.txt"
    output:
        vcf = expand("results/imputedVCFs/{{chrom}}_{i}.vcf", i=lambda wildcards: config['chunks'][wildcards.chrom])
    log:
        "logs/glimpseImpute/{chrom}.log"
    shell:
        """
        while IFS="" read -r LINE || [ -n "$LINE" ];
        do
        printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
        IRG=$(echo $LINE | cut -d" " -f3)
        ORG=$(echo $LINE | cut -d" " -f4)

        workflow/scripts/GLIMPSE/static_bins/GLIMPSE_phase_static --input {input.vcf} --reference {input.HapPanel} --map {input.genMap} --input-region ${IRG} \ 
        --output-region ${ORG} --output {output.vcf}

        bcftools index {output.vcf} 
        done < {input.chunks}
        """

 #rule glimpseLigate:
  #   input:
#
 #    output:
  #   log:
   #  params:
    # shell:
     #    """
   #      GLIMPSE_ligate --input ${LST} --output $OUT
  #       """

# rule glimpsePhase:
#     input:
#     output:
#     log:
#     params:
#     shell:
#         """
#         GLIMPSE_sample --input ${VCF} --solve --output ${OUT}
#         """


# rule glimpseConcordance:
#     input:
#     output:
#     log:
#     params:
#     shell:
#         """
#         GLIMPSE_concordance --input concordance.lst --minDP 8 --output GLIMPSE_concordance/output --minPROB 0.9999 \
#  > --bins 0.00000 0.00100 0.00200 0.00500 0.01000 0.05000 0.10000 0.20000 0.50000ex \
#  > --info_af AF_nfe --thread 4 
#         """
