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


def get_local_memory():
    '''
    Get cellranger mkref memory options.
    '''
    jobmode = config['cellranger']['jobmode']
    threads = config['cellranger']['threads']
    
    memory_value = config['cellranger']['memory_per_cpu'] * threads
    
    return memory_value


def get_runtime_options():
    '''
    Get cellranger runtime options.
    '''
    option_str = ""
    
    jobmode = config['cellranger']['jobmode']
    threads = config['cellranger']['threads']
    
    if jobmode == "local":
        local_memory = get_local_memory(),
        option_str += f"--jobmode=local --localcores={threads} --localmem={local_memory}"
    elif jobmode == "sge":
        memory_per_cpu = config['cellranger']['memory_per_cpu']
        option_str += f"--jobmode=sge --maxjobs={threads} --mempercore={memory_per_cpu}"
    else:
        raise NameError(f"Invalid job mode: {jobmode}")
    
    return option_str


def get_rule_threads():
    '''
    Get the number of threads given to run the rule.
    '''
    threads = 1
    
    jobmode = config['cellranger']['jobmode']
    if jobmode == "local":
        threads = config['cellranger']['threads']
    elif jobmode == "sge":
        threads = 1
    else:
        raise NameError(f"Invalid job mode: {jobmode}")
    
    return threads


rule cellranger_count:
    input:
        fastqs=lambda wildcards: samples["fastqs"][wildcards.sample]
    output:
        raw="{OUTDIR}/{sample}/cellranger_count/raw_feature_bc_matrix.h5",
        filtered="{OUTDIR}/{sample}/cellranger_count/filtered_feature_bc_matrix.h5",
        bam="{OUTDIR}/{sample}/cellranger_count/possorted_genome_bam.bam"
    params:
        transcriptome=config['cellranger']['transcriptome'],
        expect_cells=lambda wildcards, input: samples['expect_cells'][wildcards.sample],
        runtime_options=get_runtime_options(),
        sample_option=get_sample_option,
        threads=config['cellranger']['threads']
    envmodules:
        "bio/cellranger/3.1.0"
    threads: get_rule_threads()
    resources:
        mem_free_gb=config['cellranger']['memory_per_cpu']
    log:
        err="{OUTDIR}/logs/cellranger_count/{sample}.err",
        out="{OUTDIR}/logs/cellranger_count/{sample}.out",
        time="{OUTDIR}/logs/time/cellranger_count/{sample}"
    shell:
        """
        {DATETIME} > {log.time} &&
        rm -rf "{OUTDIR}/{wildcards.sample}/cellranger_count/" &&
        cellranger count --id={wildcards.sample} \
        --transcriptome={params.transcriptome} \
        --fastqs={input.fastqs} \
        {params.sample_option} \
        --expect-cells={params.expect_cells} \
        {params.runtime_options} \
        2> {log.err} > {log.out} &&
        mv {wildcards.sample} "{OUTDIR}/{wildcards.sample}/cellranger_count/" &&
        {DATETIME} >> {log.time}
        """
