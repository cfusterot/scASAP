import glob

FASTQDIR=samples["fastqs"]
PREFIX=samples["prefix"]
SAMPLES=samples["sample"]
read=['1', '2', '3']
lane=config["cellranger"]["lanes"]

rule fastqc:
    input:
        expand("{fqdir}/{prefix}/{prefix}_S1_L00{lane}_R{read}_001.fastq.gz", fqdir=FASTQDIR, prefix=PREFIX, read=read, lane=lane)
    output:
        expand(["{{OUTDIR}}/qc/{sample}/{prefix}_S1_L00{lane}_R{read}_001.html",
        "{{OUTDIR}}/qc/{sample}/{prefix}_S1_L00{lane}_R{read}_001_fastqc.zip"], lane=lane, prefix=PREFIX, sample=SAMPLES, read=read)
# the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: expand("--quiet --outdir '{{OUTDIR}}/qc/{sample}'", sample=SAMPLES) 
    resources:
        mem_mb=get_resource("fastqc", "mem_mb"),
        walltime=get_resource("fastqc", "walltime")
    log:
        err=expand("{{OUTDIR}}/logs/qc/{sample}_{prefix}_L00{lane}_R{read}.err", lane=lane, sample=SAMPLES, prefix=PREFIX, read=read),
        out=expand("{{OUTDIR}}/logs/qc/{sample}_{prefix}_L00{lane}_R{read}.out", lane=lane, sample=SAMPLES, prefix=PREFIX, read=read)
    threads: 
        threads=get_resource("fastqc", "threads")
    wrapper:
        "v1.23.1/bio/fastqc"
