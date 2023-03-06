import glob

rule amulet:
    input:
        finish="{OUTDIR}/{sample}/cellranger_count/cellranger.finish",
        fragments="{OUTDIR}/{sample}/cellranger_count/outs/fragments.tsv.gz",
        singlecells="{OUTDIR}/{sample}/cellranger_count/outs/singlecell.csv"
    output:
        multiplets="{OUTDIR}/{sample}/amulet/MultipletSummary.txt"
    conda:
        "../envs/amulet.yaml"
    params:
        amulet="resources/AMULET",
        autosomes=config['amulet']['autosomes'],
        blacklist=config['amulet']['blacklist']
    threads: get_resource("amulet", "threads")
    resources:
        mem_mb=get_resource("amulet", "mem_mb"),
        walltime=get_resource("amulet", "walltime")
    log:
        err="{OUTDIR}/logs/{sample}/amulet.err",
        out="{OUTDIR}/logs/{sample}/amulet.out",
    shell:
        """
        chmod +x {params.amulet}/AMULET.sh &&
        {params.amulet}/AMULET.sh {input.fragments} {input.singlecells} \
        {params.autosomes} {params.blacklist} \
        {OUTDIR}/{wildcards.sample}/amulet {params.amulet}
        """
