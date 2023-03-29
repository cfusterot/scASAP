import glob

rule mgatk:
    input:
        finish="{OUTDIR}/{sample}/cellranger_count/cellranger.finish"
    params:
        bam="{OUTDIR}/{sample}/cellranger_count/outs/possorted_bam.bam"
    output:
        ref="{OUTDIR}/{sample}/mgatk/final/chrM_refAllele.txt"
    conda:
        "../envs/mgatk.yaml"
    threads: get_resource("mgatk", "threads")
    resources:
        mem_mb=get_resource("mgatk", "mem_mb"),
        walltime=get_resource("mgatk", "walltime")
    log:
        err="{OUTDIR}/logs/{sample}/mgatk.err",
        out="{OUTDIR}/logs/{sample}/mgatk.out",
    shell:
        """
        mgatk tenx -i {params.bam} \
        -n {wildcards.sample} -o {OUTDIR}/{wildcards.sample}/mgatk/ \
        -bt CB -b {OUTDIR}/{wildcards.sample}/cellranger_count/outs/filtered_peak_bc_matrix/barcodes.tsv \ 
        2> {log.err} > {log.out} &&
        """
