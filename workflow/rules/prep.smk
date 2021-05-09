
rule referencePanel:
    input:
        vcf = "/home/sanj/ag1000g/data/ag1000g.phase2.ar1.pass.biallelic.{chrom}.vcf.gz"
    output:
        "resources/ag1000g_WestAfrica_col_{chrom}.vcf.gz"
    log:
        view = "logs/bcftoolsView/{chrom}.log",
        anno = "logs/bcftoolsAnnotate/{chrom}.log"
    params:
        samples = "resources/list_samples/WestAfrica_col_sample.list"
    shell:
        "bcftools view -Ov -S {params.samples} {input.vcf} 2> {log.view} | bcftools annotate -x INFO -Oz -o {output} 2> {log.anno}"


rule GenomeIndex:
    input:
        ref = lambda wildcards: config['ref']
    output:
        idx = touch("resources/reference/.bwa.index")
    shell:
        """
        bwa index {input.ref}
        """

rule alignBWA:
    """
    Align with bwa mem, and sorting by coordinate with samtools sort. then index with samtools.  
    """
    input:
        reads = expand("resources/reads/{{sample}}_{n}.fq.gz", n=[1,2]),
        ref = lambda wildcards: config['ref'],
        idx = "resources/reference/.bwa.index"
    output:
        bam = "resources/alignments/{sample}.bam"
    log:
        align="logs/align_bwa/{sample}.log",
        sort="logs/sort/{sample}.log",
    resources:bwa=1
    threads:8 # each sample tiny so perhaps better to run each on single thread
    params:
        tag="'@RG\\tID:{sample}\\tSM:{sample}\\tPU:nye\\tPL:nye\\tLB:{sample}_lb{sample}'"
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.reads} -R {params.tag} 2> {log.align} |
        samtools sort -@{threads} -o {output} 2> {log.sort}
        """

rule indexBams:
     input:
        "resources/alignments/{sample}.bam"
     output:
        "resources/alignments/{sample}.bam.bai"
     log:
        "logs/index_bams/{sample}.log"
     shell:
        "samtools index {input} {output} 2> {log}"


rule lowCovGenotypeLikelihoods:
    """
    Get pileup of reads at target loci and pipe output to bcftoolsCall
    """
    input:
        bam = "resources/alignments/{sample}.bam",
        index = "resources/alignments/{sample}.bam.bai",
        ref = lambda wildcards: config['ref'],
    output:
        calls = "results/vcfs/{sample}.calls.vcf"
    log:
        mpileup = "logs/mpileup/{sample}.log",
        call = "logs/bcftools_call/{sample}.log"
    params:
        regions = "haplptype_ref_sites_vcf", #*************
        depth = 2000
    shell:
        """
        bcftools mpileup -Ou -f {params.ref} -I -E -R {params.regions} -a 'FORMAT/DP' -T {params.regions} {input.bam} 2> {log.mpileup} |
        bcftools call -Aim -C alleles -T {params.regions} -Ou 2> {log.call} | bcftools sort -Oz -o {output.calls} 2> {log.call}
        """













# rule mergeVCFs:




