



rule GenomeIndexCOE:
    input:
        ref = lambda wildcards: config['coeref']
    output:
        idx = touch("results/.bwa.index.coe")
    shell:
        """
        bwa index {input.ref}
        """


rule alignBWA_coe:
    """
    Align with bwa mem, and sorting by coordinate with samtools sort. then index with samtools.  
    """
    input:
        reads = expand("resources/reads/DW-{{sample}}_R{n}_001.fastq.gz", n=[1,2]),
        ref = lambda wildcards: config['coeref'],
        idx = "results/.bwa.index.coe"
    output:
        bam = "results/alignments/coe/{sample}.bam"
    log:
        align = "logs/align_bwa/{sample}_coe.log",
        sort = "logs/sort/{sample}_coe.log",
    resources: bwa = 1 
    threads:8 # each sample tiny so perhaps better to run each on single thread
    params:
        tag="'@RG\\tID:{sample}\\tSM:{sample}\\tPU:nye\\tPL:nye\\tLB:{sample}_lb{sample}'"
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.reads} -R {params.tag} 2> {log.align} |
        samtools sort -@{threads} -o {output} 2> {log.sort}
        """

rule indexBams_coe:
     input:
        "results/alignments/coe/{sample}.bam"
     output:
        "results/alignments/coe/{sample}.bam.bai"
     log:
        "logs/index_bams_coe/{sample}.log"
     shell:
        "samtools index {input} {output} 2> {log}"



rule lowCovGenotypeLikelihoods_coe:
    """
    Get pileup of reads at target loci and pipe output to bcftoolsCall
    """
    input:
        bam = "results/alignments/coe/{sample}.bam",
        index = "results/alignments/coe/{sample}.bam.bai",
        vcf = "resources/ag1000g_WestAfrica_col_2L.sites.vcf.gz",
        tsv = "resources/ag1000g_WestAfrica_col_2L.sites.tsv.gz",
        ref = config['coeref'],
    output:
        calls = "results/coevcfs/{sample}.calls.28mb.vcf.gz"
    log:
        mpileup = "logs/mpileup_coe/{sample}.log",
        call = "logs/bcftools_call_coe/{sample}.log"
    shell:
        """
        bcftools mpileup -Ou -f {input.ref} -I -E -a 'FORMAT/DP' -T {input.vcf} -r 2L {input.bam} 2> {log.mpileup} |
        bcftools call -Aim -C alleles -T {input.tsv} -Ou 2> {log.call} | bcftools sort -Oz -o {output.calls} 2> {log.call}
        """