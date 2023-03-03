import glob

rule mgatk:
    input:
        finish="{OUTDIR}/{sample}/cellranger_count/cellranger.finish",
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
        err="{OUTDIR}/logs/mgatk/{sample}.err",
        out="{OUTDIR}/logs/mgatk/{sample}.out",
        time="{OUTDIR}/logs/time/mgatk/{sample}"
    shell:
        """
        {DATETIME} > {log.time} &&
        mgatk tenx -i {input.bam} \
        -n {wildcards.sample} -o {OUTDIR}/{wildcards.sample}/mgatk/ \
        -bt CB -b {OUTDIR}/{wildcards.sample}/cellranger_count/outs/filtered_peak_bc_matrix/barcodes.tsv \ 
        2> {log.err} > {log.out} &&
        {DATETIME} >> {log.time}
        """
