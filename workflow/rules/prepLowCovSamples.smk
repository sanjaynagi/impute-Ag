rule genomeIndex:
    """
    Index the reference genome with BWA
    """
    input:
        ref = config['ref']
    output:
        idx = touch("resources/reference/.bwa.index")
    log:
        "logs/align_bwa/index.log"
    shell:
        """
        bwa index {input.ref} 2> {log}
        """

rule alignBWA:
    """
    Align with bwa mem, and sorting by coordinate with samtools sort. then index with samtools.  
    """
    input:
        reads = expand("resources/reads/{{sample}}_{n}.fq.gz", n=[1,2]),
        ref = config['ref'],
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
    Rule that downsamples a bam file to a required number of reads, if less reads then dont try and downsample
    """
    input:
        bam = "results/alignments/{sample}.bam"
    output:
        reducedBam = "results/alignments/downSampled/{sample}.bam",
#	index = "results/alignments/downSampled/{sample}.bam.bai"
    log:
        "logs/subSampleBams/{sample}.log"
    params:
        reads = config['subsample']['reads']
    shell:
        """
        FACTOR=$(samtools idxstats {input.bam} | cut -f3 | awk -v COUNT={params.reads} 'BEGIN {{total=0}} {{total += $1}} END {{print COUNT/total}}')

        if [[ $FACTOR > 1 ]]
        then 
        echo '[ERROR]: Requested number of reads exceeds total read count in' {input.bam}, not downsampling
        cp {input.bam} {output.reducedBam} && \
        samtools index {output.reducedBam} {output.reducedBam}.bai 2>> {log}
        fi
    
        if [[ $FACTOR < 1 ]]
        then
        sambamba view -s $FACTOR -f bam -l 5 {input.bam} > {output.reducedBam} 2> {log} && \
        samtools index {output.reducedBam} {output.reducedBam}.bai 2>> {log}
        fi
        """



rule lowCovGenotypeLikelihoods:
    """
    Get pileup of reads at target loci and pipe output to bcftoolsCall
    """
    input:
        bam = whichBams(),
        index = lambda wildcards: whichBams() + ".bai",
        vcf = "resources/ag3.{chrom}.sites.vcf.gz",
        csi = "resources/ag3.{chrom}.sites.vcf.gz.csi",
        tsv = "resources/ag3.{chrom}.sites.tsv.gz",
        tbi = "resources/ag3.{chrom}.sites.tsv.gz.tbi",
        ref = config['ref'],
    output:
        calls = "results/vcfs/{sample}.calls.{chrom}.vcf.gz"
    log:
        mpileup = "logs/mpileup/{sample}_{chrom}.log",
        call = "logs/bcftoolsCall/{sample}_{chrom}.log"
    shell:
        """
        bcftools mpileup -Ou -f {input.ref} -I -E -a 'FORMAT/DP' -T {input.vcf} -r {wildcards.chrom} {input.bam} 2> {log.mpileup} |
        bcftools call -Aim -C alleles -T {input.tsv} -Ou 2> {log.call} | bcftools sort -Oz -o {output.calls} 2>> {log.call}
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
        expand("results/vcfs/{sample}.calls.{{chrom}}.vcf.gz", sample=samples),
        expand("results/vcfs/{sample}.calls.{{chrom}}.vcf.gz.csi", sample=samples)
     output:
        "results/vcfs/merged_calls.{chrom}.vcf.gz"
     log:
        "logs/mergeVCFs/{chrom}.log"
     shell:
        """
        find results/{wildcards.dataset}VCFs/ -name *calls.{wildcards.chrom}.vcf.gz | sort > results/{wildcards.dataset}VCFs/sampleVCF.{wildcards.chrom}.list 2> {log}
        bcftools merge -m none -r {wildcards.chrom} -Oz -o {output} -l results/{wildcards.dataset}VCFs/sampleVCF.{wildcards.chrom}.list 2>> {log}
        """

rule indexMergedVCFs:
     input:
        "results/vcfs/merged_calls.{chrom}.vcf.gz"
     output:
        "results/vcfs/merged_calls.{chrom}.vcf.gz.csi"
     log:
        "logs/indexMergedVCFs/{chrom}.log"
     shell:
        "bcftools index {input} 2> {log}"
