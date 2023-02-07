import glob

rule fastqc:
    input:
        lambda wc: "samples["fastqs"]/{sample}_R1_001.fastq.gz",
        lambda wc: "samples["fastqs"]/{sample}_R2_001.fastq.gz"
        #r1=lambda wildcard: glob.glob(FASTQDIR + [wildcards.sample] + '_R1.fastq.gz'),
        #r2=lambda wildcard: glob.glob(FASTQDIR + [wildcards.sample] + '_R2.fastq.gz')
    output:
        html="{OUTDIR}/qc/{sample}/{sample}.html",
        zip="{OUTDIR}/qc/{sample}/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "--quiet --outdir {OUTDIR}/qc/"
    resources:
        mem_mb=get_resource("default", "mem_mb"),
        walltime=get_resource("default", "walltime")
    log:
        err="{OUTDIR}/logs/qc/{sample}.err",
        out="{OUTDIR}/logs/qc/{sample}.out"
    threads: 1
    wrapper:
        "v1.23.1/bio/fastqc"
