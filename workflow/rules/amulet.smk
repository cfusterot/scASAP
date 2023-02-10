import glob

rule amulet:
    input:
        
    output:
        multiplets="{OUTDIR}/amulet/{sample}/MultipletSummary.txt"
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
