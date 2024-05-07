import glob

rule cnv:
    input:
        finish="{}/{{sample}}/cellranger_count/cellranger.finish".format(OUTDIR)
    output:
        "{}/{{sample}}/cnv/cnv_calls.rds".format(OUTDIR)
    conda:
        "../envs/epianeufinder.yaml"
    params:
        fragments="{}/{{sample}}/cellranger_count/fragments.tsv.gz".format(OUTDIR),
        blacklist=config['amulet']['blacklist'],
        windowSize = config['cnv']['windowSize'],
        genome = config['cnv']['genome'],
        exclude = config['cnv']['exclude'],
        output_dir="{}/{{sample}}/cnv".format(OUTDIR),
    threads: get_resource("cnv", "threads")
    resources:
        mem_mb=get_resource("cnv", "mem_mb"),
        walltime=get_resource("cnv", "walltime")
    log:
        err="{}/{{sample}}/amulet.err".format(LOGDIR),
        out="{}/{{sample}}/amulet.out".format(LOGDIR)
    script: 
        "../scripts/epianeufinder.R"
