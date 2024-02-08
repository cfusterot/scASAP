import glob

rule amulet:
    input:
        finish="{}/{{sample}}/cellranger_count/cellranger.finish".format(OUTDIR)
    output:
        multiplets="{}/{{sample}}/amulet/MultipletSummary.txt".format(OUTDIR)
    conda:
        "../envs/amulet.yaml"
    params:
        fragments="{}/{{sample}}/cellranger_count/fragments.tsv.gz".format(OUTDIR),
        singlecells="{}/{{sample}}/cellranger_count/singlecell.csv".format(OUTDIR),
        amulet="resources/AMULET",
        autosomes=config['amulet']['autosomes'],
        blacklist=config['amulet']['blacklist']
    threads: get_resource("amulet", "threads")
    resources:
        mem_mb=get_resource("amulet", "mem_mb"),
        walltime=get_resource("amulet", "walltime")
    log:
        err="{}/{{sample}}/amulet.err".format(LOGDIR),
        out="{}/{{sample}}/amulet.out".format(LOGDIR)
    shell:
        """
        chmod +x {params.amulet}/AMULET.sh &&
        {params.amulet}/AMULET.sh {params.fragments} {params.singlecells} \
        {params.autosomes} {params.blacklist} \
        {OUTDIR}/{wildcards.sample}/amulet {params.amulet}
        """
