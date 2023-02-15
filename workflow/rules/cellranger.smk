import glob
import os

def get_sample_option(wc):
    option_str=""
    prefix=units["sample"][wc.sample] 
    option_str += f"--sample={prefix}"
    return option_str

def get_fq_path(wc):
    p=units["fq1"][wc.sample]
    return p

rule cellranger_count:
    input:
        fq=os.path.dirname(units['fq1'])
        #fq=os.path.dirname(get_fq_path)
        #fq= wildcards: os.path.dirname(units["fq1"][wildcards.sample])
        #fq2=lambda wildcards: units["fq2"][wildcards.sample],
        #fq3=lambda wildcards: units["fq3"][wildcards.sample]
    output:
        raw="{OUTDIR}/cellranger_count/{sample}/raw_feature_bc_matrix.h5",
        filtered="{OUTDIR}/cellranger_count/{sample}/filtered_feature_bc_matrix.h5",
        bam="{OUTDIR}/cellranger_count/{sample}/possorted_genome_bam.bam"
    params:
        reference=config['cellranger']['reference'],
        sample_option=get_sample_option
    envmodules:
        config['envmodules']['cellranger']
    threads: get_resource("cellranger", "threads")
    resources:
        mem_mb=get_resource("cellranger", "mem_mb"),
        walltime=get_resource("cellranger", "walltime")
    log:
        err="{OUTDIR}/logs/cellranger_count/{sample}.err",
        out="{OUTDIR}/logs/cellranger_count/{sample}.out",
        time="{OUTDIR}/logs/time/cellranger_count/{sample}"
    shell:
        """
        {DATETIME} > {log.time} &&
        cellranger-atac count --id={wildcards.sample} \
        --reference={params.reference} \
        --fastqs={input.fq} \
        {params.sample_option} \
        2> {log.err} > {log.out} &&
        mv {wildcards.sample} "{OUTDIR}/cellranger_count/" &&
        {DATETIME} >> {log.time}
        """
