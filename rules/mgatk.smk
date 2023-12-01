import glob

rule mgatk:
    input:
        finish="{}/{{sample}}/cellranger_count/cellranger.finish".format(OUTDIR)
    params:
        bam="{}/{{sample}}/cellranger_count/outs/possorted_bam.bam".format(OUTDIR)
    output:
        ref="{}/{{sample}}/mgatk/final/{{sample}}.variant_stats.tsv.gz".format(OUTDIR)
    conda:
        "../envs/mgatk.yaml"
    threads: get_resource("mgatk", "threads")
    resources:
        mem_mb=get_resource("mgatk", "mem_mb"),
        walltime=get_resource("mgatk", "walltime")
    log:
        err="{}/{{sample}}/mgatk.err".format(LOGDIR),
        out="{}/{{sample}}/mgatk.out".format(LOGDIR)
    shell:
        """
        # delete folder if exists
        rm -f -R {OUTDIR}/{wildcards.sample}/mgatk
        # common error resolved by those two export commands
        export LC_ALL=C.UTF-8
        export LANG=C.UTF-8
        # run mgatk command
        mgatk tenx -i {params.bam} -n {wildcards.sample} -o {OUTDIR}/{wildcards.sample}/mgatk/ -bt CB -b {OUTDIR}/{wildcards.sample}/cellranger_count/outs/filtered_peak_bc_matrix/barcodes.tsv 
        """
