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


rule subSampleBam:
    """
    Rule that downsamples a bam file to a required number of reads
    """
    input:
        bam = "results/alignments/{sample}.bam"
    output:
        reducedBam = "results/alignments/downSampled/{sample}.bam"
    log:
        "logs/subSampleBams/{sample}.log"
    params:
        reads = 5560000
    shell:
        """
        FACTOR=$(samtools idxstats {input.bam} | cut -f3 | awk -v COUNT=$2 'BEGIN {total=0} {total += $1} END {print COUNT/total}')

        if [[ $FACTOR > 1 ]]
        then 
        echo '[ERROR]: Requested number of reads exceeds total read count in' $1 '-- exiting' && exit 1
        fi

        sambamba view -s $FACTOR -f bam -l 5 $1 > {output.reducedBam}
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

