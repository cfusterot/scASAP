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
        amulet_exec="resources/AMULET/AMULET.sh",
        amulet_folder="resources/AMULET",
        autosomes=config['amulet']['autosomes'],
        blacklist=config['amulet']['blacklist']
    threads: get_resource("amulet", "threads")
    resources:
        mem_mb=get_resource("amulet", "mem_mb"),
        walltime=get_resource("amulet", "walltime")
    log:
        err="{OUTDIR}/logs/amulet/{sample}.err",
        out="{OUTDIR}/logs/amulet/{sample}.out",
    shell:
        """
        {params.amulet_exec} {input.fragments} {input.singlecells} \ 
        {params.autosomes} {params.blacklist} {OUTDIR}/{wildcards.sample}/amulet \
        {params.amulet_folder}
        """
