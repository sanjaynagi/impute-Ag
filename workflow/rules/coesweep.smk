



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
        bam = "results/COEalignments/{sample}.bam"
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
        samtools sort -@{threads} 2> {log.sort} | samtools view -F 4 -o {output.bam} 2>> {log.sort}
        """

rule indexBams_coe:
     input:
        "results/COEalignments/{sample}.bam"
     output:
        "results/COEalignments/{sample}.bam.bai"
     log:
        "logs/index_COEbams/{sample}.log"
     shell:
        "samtools index {input} {output} 2> {log}"





rule mpileupIR:
    """
    Get allele count tables of variants of choice (specified in config file ("IRmutations.tsv"))
    """
    input:
        bam="results/alignments/{sample}.bam",
        index="results/alignments/{sample}.bam.bai",
    output:
        "results/allele_balance/counts/{sample}_{mut}_allele_counts.tsv",
    conda:
        "../envs/variants.yaml"
    priority: 10
    log:
        "logs/mpileupIR/{sample}_{mut}.log",
    params:
        region=lambda wildcards: mutationData[
            mutationData.Name == wildcards.mut
        ].Location.tolist(),
        ref=config["ref"]["genome"],
        basedir=workflow.basedir,
    shell:
        """
        samtools mpileup {input.bam} -r {params.region} -f {params.ref} 2> {log} | 
        python2 {params.basedir}/scripts/BaseParser.py > {output} 2>> {log}
        """


rule AlleleBalanceIR:
    """
    R script to take allele count tables from mpileupIR rule and output .xlsx report for all mutations of interest
    """
    input:
        counts=expand(
            "results/allele_balance/counts/{sample}_{mut}_allele_counts.tsv",
            sample=coesamples,
            mut=mutationData.Name,
        ),
        samples=coesamples,
        mutations="resources/coemutations.tsv",
    output:
        expand(
            "results/allele_balance/csvs/{mut}_allele_balance.csv",
            mut=mutationData.Name,
        ),
        allele_balance="results/allele_balance/allele_balance.xlsx",
        mean_allele_balance="results/allele_balance/mean_allele_balance.xlsx",
    conda:
        "../envs/diffexp.yaml"
    priority: 10
    log:
        "logs/AlleleBalanceIR.log",
    script:
        "../scripts/MutAlleleBalance.R"





















rule lowCovGenotypeLikelihoods_coe:
    """
    Get pileup of reads at target loci and pipe output to bcftoolsCall
    """
    input:
        bam = "results/COEalignments/{sample}.bam",
        index = "results/COEalignments/{sample}.bam.bai",
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
        bcftools mpileup -Ou -f {input.ref} -I -E -a 'FORMAT/DP' -T {input.vcf} {input.bam} 2> {log.mpileup} |
        bcftools call -Aim -C alleles -T {input.tsv} -Ou 2> {log.call} | bcftools sort -Oz -o {output.calls} 2> {log.call}
        """
