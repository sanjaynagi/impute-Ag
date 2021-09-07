
#rule subsetReferencePanel:
#    """
 #   Subset the reference panel to a specific set of samples. Only necessary to save space - runtime should be lower with larger Hap Panel
 #   """
  #  input:
  #      vcf = "/home/sanj/ag1000g/data/ag3.ar1.pass/ag3.ar1.pass.{chrom}.vcf.gz"
#    output:
#        "resources/ag1000g_WestAfrica_col_{chrom}.vcf.gz"
#    log:
#        view = "logs/bcftoolsView/{chrom}.log",
#        anno = "logs/bcftoolsAnnotate/{chrom}.log"
#    params:
#        samples = "resources/list_samples/WestAfrica_col_sample.list"
#    shell:
#        "bcftools view -Ov -S {params.samples} {input.vcf} 2> {log.view} | bcftools annotate -x INFO -Oz -o {output} 2> {log.anno}"


rule indexHapPanel:
     input:
        vcf=config['ag3']['vcf']
     output:
        csi=config['ag3']['vcf'] + ".csi"
     log:
        "logs/indexHapPanel.{chrom}.log"
     shell:
        "bcftools index {input} 2> {log}"


rule validationHapPanel:
    input:
        vcf = config['ag3']['vcf']
    output:
        vcf = "resources/ag3.validation.vcf"
    log:
        "logs/validationHap.log"
    params:
        validationSamples = config['validate']['samplesRemove']
    threads: 8
    shell:
        """
        bcftools norm -m -any {input.vcf} -Ou --threads {threads} | 
        bcftools view -m 2 -m 2 -v snps -s ^{params.validationSamples} --threads {threads} -Ob -o {output.vcf}
        bcftools index -f {output.vcf} --threads {threads}
        """

rule validationSampleVCF:
    input:
        vcf=config['ag3']['vcf']
    output:
        vcf= "results/ag3.validation.vcf.gz"
    log:
        "logs/validationVCF.log"
    params:
        validationSamples = config['validate']['samplesRemove']
    shell:
        """
        bcftools view -S {params.validationSamples} {input.vcf} -Oz -o {output.vcf} 2> {log}
        """

rule extractSites:
    """
    Extract variable sites from the Haplotype reference panel
    """
    input:
        "resources/ag3.validation.vcf" if config['validate']['activate'] else config['ag3']['vcf']
    output:
        "resources/ag3.{chrom}.sites.vcf.gz"
    log:
        "logs/extractSites.{chrom}.log"
    shell:
        """
        bcftools view -G -m 2 -M 2 -v snps {input} -Oz -o {output} 2> {log}
        """

rule indexSites:
     input:
        vcf="resources/ag3.{chrom}.sites.vcf.gz"
     output:
        csi="resources/ag3.{chrom}.sites.vcf.gz.csi"
     log:
        "logs/indexSites.{chrom}.log"
     shell:
        "bcftools index {input} 2> {log}"

rule sitesTable:
    """
    Convert sites VCF to a TSV
    """
    input:
        vcf="resources/ag3.{chrom}.sites.vcf.gz",
        csi="resources/ag3.{chrom}.sites.vcf.gz.csi"
    output:
        tsv="resources/ag3.{chrom}.sites.tsv.gz"
    log:
        "logs/sitesTable.{chrom}.log"
    shell:
        """
        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {input.vcf} | bgzip -c > {output.tsv} 2> {log}
        """

rule tabixTable:
    input:
        tsv="resources/ag3.{chrom}.sites.tsv.gz"
    output:
        tbi="resources/ag3.{chrom}.sites.tsv.gz.tbi",
    log:
        "logs/tabixTable.{chrom}.log"
    shell:
        """
        tabix -s1 -b2 -e2 {input.tsv} 2> {log}
        """
