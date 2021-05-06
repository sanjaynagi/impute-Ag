
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



# # we need to get genotype likelihoods
# rule lowCovSequencing:









# rule GenotypeLikelihoods:




# rule mergeVCFs:




