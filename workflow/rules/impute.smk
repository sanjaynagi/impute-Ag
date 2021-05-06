# An example collection of Snakemake rules imported in the main Snakefile.

# rule glimpseChunk:
#     input:
#     output:
#     log:
#     params:
#     shell:
#         """
#         GLIMPSE_chunk --input reference_panel/1000GP.chr22.noNA12878.sites.vcf.gz --region chr22 --window-size 2000000 \
#         --buffer-size 200000 --output chunks.chr22.txt 
#         """

# rule glimpseImpute:
#     input:
#     output:
#     log:
#     params:
#     shell:
#         """
#         GLIMPSE_phase --input ${VCF} --reference ${REF} --map ${MAP} --input-region ${IRG} \ 
#         --output-region ${ORG} --output ${OUT}
#         """

# rule glimpseLigate:
#     input:
#     output:
#     log:
#     params:
#     shell:
#         """
#         GLIMPSE_ligate --input ${LST} --output $OUT
#         """

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
