
rule subsetReferencePanel:
    """
    Subset the reference panel to a specific set of samples. Only necessary to save space - runtime should be lower with larger Hap Panel
    """
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
    """
    Index the reference genome with BWA
    """
    input:
        ref = config['ref']
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
        bam = "results/alignments/{sample}.bam"
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
        "results/alignments/{sample}.bam"
     output:
        "results/alignments/{sample}.bam.bai"
     log:
        "logs/indexBams/{sample}.log"
     shell:
        "samtools index {input} {output} 2> {log}"

rule extractSites:
    """
    Extract variable sites from the Haplotype reference panel
    """
    input:
        "resources/ag1000g_phase2.{chrom}.vcf.gz"
    output:
        "resources/ag1000g_phase2.{chrom}.sites.vcf.gz"
    log:
        "logs/extractSites.{chrom}.log"
    shell:
        """
        bcftools view -G -m 2 -M 2 -v snps {input} -0z -o {output} 2> {log}
        """

rule indexSites:
     input:
        vcf="resources/ag1000g_phase2.{chrom}.sites.vcf.gz"
     output:
        csi="resources/ag1000g_phase2.{chrom}.sites.vcf.gz.csi"
     log:
        "logs/indexSites.{chrom}.log"
     shell:
        "bcftools index {input} 2> {log}"

rule sitesTable:
    """
    Convert sites VCF to a TSV
    """
    input:
        vcf="resources/ag1000g_phase2.{chrom}.sites.vcf.gz"
    output:
        tsv="resources/ag1000g_phase2.{chrom}.sites.tsv.gz"
    log:
        "logs/sitesTable.{chrom}.log"
    shell:
        """
        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {input.vcf} | bgzip -c > {output.tsv} 2> {log}
        """

rule tabixTable:
    input:
        tsv="resources/ag1000g_phase2.{chrom}.sites.tsv.gz"
    output:
        tbi="resources/ag1000g_phase2.{chrom}.sites.tsv.gz.tbi",
    log:
        "logs/tabixTable.{chrom}.log"
    shell:
        """
        tabix -s1 -b2 -e2 {output.tsv} 2> {log}
        """

rule lowCovGenotypeLikelihoods:
    """
    Get pileup of reads at target loci and pipe output to bcftoolsCall
    """
    input:
        bam = "results/alignments/{sample}.bam",
        index = "results/alignments/{sample}.bam.bai",
        vcf = "resources/ag1000g_phase2.{chrom}.sites.vcf.gz",
        tsv = "resources/ag1000g_phase2.{chrom}.sites.tsv.gz",
        ref = config['ref'],
    output:
        calls = "results/vcfs/{sample}.calls.{chrom}.vcf.gz"
    log:
        mpileup = "logs/mpileup/{sample}_{chrom}.log",
        call = "logs/bcftoolsCall/{sample}_{chrom}.log"
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

