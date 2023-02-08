import glob


FASTQDIR=samples["fastqs"]
PREFIX=samples["prefix"]
SAMPLES=samples["sample"]

rule fastqc:
    input:
        expand("{fqdir}/{prefix}/{prefix}_S1_L002_I1_001.fastq.gz", fqdir=FASTQDIR, prefix=PREFIX)      
    output:
        html=expand("{{OUTDIR}}/qc/{sample}/{prefix}.html", sample=SAMPLES, prefix=PREFIX),
        zip=expand("{{OUTDIR}}/qc/{sample}/{prefix}_fastqc.zip", sample=SAMPLES, prefix=PREFIX) # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "--quiet --outdir {OUTDIR}/qc/"
    resources:
        mem_mb=get_resource("default", "mem_mb"),
        walltime=get_resource("default", "walltime")
    log:
        err=expand("{{OUTDIR}}/logs/qc/{sample}-{prefix}.err", sample=SAMPLES, prefix=PREFIX),
        out=expand("{{OUTDIR}}/logs/qc/{sample}-{prefix}.out", sample=SAMPLES, prefix=PREFIX)
    threads: 1
    wrapper:
        "v1.23.1/bio/fastqc"
