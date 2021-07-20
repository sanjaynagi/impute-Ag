

def whichBams():
    if config['subsample']['activate']:
        bam = "results/alignments/downSampled/{sample}.bam"
    else:
        bam = "results/alignments/{sample}.bam"

    return bam
