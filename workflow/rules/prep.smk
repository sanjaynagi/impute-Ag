
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




"""
extract variable sites from hap panel 
bcftools view -G -m 2 -M 2 -v snps ag1000g_WestAfrica_col_X.vcf.gz -Oz -o ag1000g_WestAfrica_col_X.sites.vcf.gz
bcftools index ag1000g_WestAfrica_col_2L.sites.vcf.gz

bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' ag1000g_WestAfrica_col_X.sites.vcf.gz | bgzip -c > ag1000g_WestAfrica_col_X.sites.tsv.gz
tabix -s1 -b2 -e2 ag1000g_WestAfrica_col_X.sites.tsv.gz

"""







rule lowCovGenotypeLikelihoods:
    """
    Get pileup of reads at target loci and pipe output to bcftoolsCall
    """
    input:
        bam = "resources/alignments/{sample}.bam",
        index = "resources/alignments/{sample}.bam.bai",
        vcf = "resources/ag1000g_WestAfrica_col_{chrom}.sites.vcf.gz",
        tsv = "resources/ag1000g_WestAfrica_col_{chrom}.sites.tsv.gz",
        ref = lambda wildcards: config['ref'],
    output:
        calls = "results/vcfs/{sample}.calls.{chrom}.vcf.gz"
    log:
        mpileup = "logs/mpileup/{sample}_{chrom}.log",
        call = "logs/bcftools_call/{sample}_{chrom}.log"
    shell:
        """
        bcftools mpileup -Ou -f {input.ref} -I -E -a 'FORMAT/DP' -T {input.vcf} -r {wildcards.chrom} {input.bam} 2> {log.mpileup} |
        bcftools call -Aim -C alleles -T {input.tsv} -Ou 2> {log.call} | bcftools sort -Oz -o {output.calls} 2> {log.call}
        """



rule indexVCFs:
     input:
        "results/vcfs/{sample}.calls.{chrom}.vcf.gz"
     output:
        "results/vcfs/{sample}.calls.{chrom}.vcf.gz.csi"
     log:
        "logs/indexVCFs/{sample}_{chrom}.log"
     shell:
        "bcftools index {input} 2> {log}"




rule mergeVCFs:
     input:
        expand("results/vcfs/{sample}.calls.{{chrom}}.vcf", sample=samples)
     output:
        "results/vcfs/merged_calls.{chrom}.vcf.gz"
     log:
        "logs/mergeVCFs/{chrom}.log"
     params:
        list = "resources/vcf.{chrom}.list"
     shell:
        "bcftools merge -m none -r {wildcards.chrom} -Oz -o {output} -l {params.list} 2> {log}"









# rule mergeVCFs:




