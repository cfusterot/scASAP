import glob

rule amulet:
    input:
        fragments="{OUTDIR}/{sample}/cellranger_count/{sample}/outs/fragments.tsv.gz",
        singlecells="{OUTDIR}/{sample}/cellranger_count/{sample}/outs/singlecell.csv"
    output:
        multiplets="{OUTDIR}/{sample}/amulet/MultipletSummary.txt"
    params:
        amulet_exec="/resources/AMULET/AMULET.sh",
        amulet_folder="/resources/AMULET",
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
        {params.autosomes} {params.blacklist} "{OUTDIR}/{sample}/amulet" \
        {params.amulet_folder}
        """
