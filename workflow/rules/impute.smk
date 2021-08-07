# An example collection of Snakemake rules imported in the main Snakefile.

rule glimpseChunk:
    """
    Divide the genome into chunks
    """
    input:
        sites = "resources/ag3.{chrom}.sites.vcf.gz",
    output:
        chunks = "resources/chunks.{chrom}.txt"
    log:
        "logs/glimpseChunk/{chrom}.log"
    shell:
        """
        workflow/scripts/GLIMPSE/static_bins/GLIMPSE_chunk_static --input {input.sites} --region {wildcards.chrom} --window-size 2000000 \
        --buffer-size 200000 --output {output.chunks} 2> {log}
        """

rule glimpseImputeLigate:
    """
    Run imputation and phasing algorithm
    """
    input:
        vcf = "results/{dataset}VCFs`/merged_calls.{chrom}.vcf.gz",
        index = "results/{dataset}VCFs`/merged_calls.{chrom}.vcf.gz.csi",
        HapPanel = config['ag3']['vcf'],
        genMap = "resources/geneticMaps/{chrom}.gmap",
        chunks = "resources/chunks.{chrom}.txt"
    output:
        vcf = "results/{dataset}VCFs/imputed.{chrom}.vcf.gz"
    log:
        "logs/glimpseImputeLigate/{chrom}_{dataset}.log"
    threads: 24
    shell:
        """
        workflow/scripts/glimpsePhaseLigate.sh {wildcards.chrom} {input.vcf} {input.HapPanel} {input.genMap} {input.chunks} {threads} {wildcards.dataset} 2> {log}
        """

rule glimpsePhase:
    input:
        vcf = "results/{dataset}VCFs/imputed.{dataset}.{chrom}.vcf.gz"
    output:
        phasedVCF = "results/{dataset}VCFs`/phased.{chrom}.vcf.gz"
    log:
        "logs/glimpseHaplotypes.{chrom}_{dataset}.log"
    threads: 4
    shell:
        """
        workflow/scripts/GLIMPSE/static_bins/GLIMPSE_sample_static --input {input.vcf} --solve --thread {threads} --output {output.phasedVCF} 2> {log}
        """


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
