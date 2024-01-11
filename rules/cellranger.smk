import glob
import os
from pathlib import Path

def get_sample_option(wc):
    option_str=""
    prefix=units["sample"][wc.sample] 
    option_str += f"--sample={prefix}"
    return option_str

def get_fq_path(wc):
    filepath = units["fq1"][wc.sample][0]
    dirname, basename = os.path.split(filepath)
    return dirname

rule cellranger_count:
    input:
        fq=lambda wc: units["fqs"][wc.sample]
    output:
        finish="{}/{{sample}}/cellranger.finish".format(LOGDIR)
    params:    
        reference=config['cellranger']['reference']
    envmodules:
        config['envmodules']['cellranger']
    threads: get_resource("cellranger", "threads")
    resources:
        mem_mb=get_resource("cellranger", "mem_mb"),
        walltime=get_resource("cellranger", "walltime")
    log:
        err="{}/{{sample}}/cellranger.err".format(LOGDIR),
        out="{}/{{sample}}/cellranger.out".format(LOGDIR),
    shell:
        """
<<<<<<< HEAD
        {DATETIME} > {log.time} &&
        cellranger-atac count --id={wildcards.sample} \
        --reference={params.reference} \
        --fastqs={input.fq} 2> {log.err} > {log.out}
        cp -R {wildcards.sample}/* {OUTDIR}/{wildcards.sample}/cellranger_count 
        rm -R {wildcards.sample}
        touch {OUTDIR}/{wildcards.sample}/cellranger_count/cellranger.finish
        {DATETIME} >> {log.time} 
=======
        cellranger-atac count --id={wildcards.sample} --reference={params.reference} --fastqs={input.fq} 2> {log.err} > {log.out}
        {DATETIME} >> {output.finish} 
>>>>>>> 43c662ef70b5fb76c8016abe59e047268bc5a419
        """

rule cellranger_mv:
    input:
        log="{}/{{sample}}/cellranger.finish".format(LOGDIR)
    output:
        finish="{}/{{sample}}/cellranger_count/cellranger.finish".format(OUTDIR)
    threads: get_resource("default", "threads")
    resources:
        mem_mb=get_resource("default", "mem_mb"),
        walltime=get_resource("default", "walltime")
    log:
        err="{}/{{sample}}/cellranger_mv.err".format(LOGDIR),
        out="{}/{{sample}}/cellranger_mv.out".format(LOGDIR)
    shell:
       """
       rsync -av {wildcards.sample}/* {OUTDIR}/{wildcards.sample}/cellranger_count
       rm -rf {wildcards.sample}
       touch {OUTDIR}/{wildcards.sample}/cellranger_count/cellranger.finish
       """
    
