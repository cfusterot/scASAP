import glob


def get_sample_option(wildcards):
    '''
    Get cellranger sample option.
    '''
    
    option_str = ""
    
    sample_prefix = samples['prefix'][wildcards.sample]

    if sample_prefix != ".":
        option_str += f"--sample={sample_prefix}"
        
    return option_str


rule cellranger_count:
    input:
        fastqs=lambda wildcards: samples["fastqs"][wildcards.sample]
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
        --fastqs={input.fastqs} \
        {params.sample_option} \
        --localcores=12 --localmem=40 \
        2> {log.err} > {log.out} &&
        mv {wildcards.sample} "{OUTDIR}/cellranger_count/" &&
        {DATETIME} >> {log.time}
        """
