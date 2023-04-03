import glob

rule mutcalling:
    input:
        singlecell = "{}/{{sample}}/cellranger_count/outs/singlecell.csv".format(OUTDIR)
    output:
        mut_output=dynamic("{OUTDIR}/{sample}/gotcha/results/Split/Filtered/_rslurm_{letter}/results_0.RDS")
    params:
        input_dir = samples["gotcha"][wc.sample],
        output_dir = "{OUTDIR}/{sample}/gotcha/",
        singlecell = "{}/{{sample}}/cellranger_count/outs/singlecell.csv".format(OUTDIR),
        whitelist=config['gotcha']['whitelist'],
        primer=config['gotcha']['primer']['primer_sequence'],
        reverse_complement = config['gotcha']['primer']['reverse_complement'],
        primed_max_mismatch = config['gotcha']['primer']['max_mismatch'],
        read = config['gotcha']['mutation']['read'],
        wt_sequence =  config['gotcha']['mutation'['wt_sequence'],
        mut_sequence = config['gotcha']['mutation']['mut_sequence'],
        mut_start = config['gotcha']['mutation'['start'],
        mut_end = config['gotcha']['mutation'['end'],
    threads: get_resource("gotcha", "threads")
    resources:
        mem_mb=get_resource("gotcha", "mem_mb"),
        walltime=get_resource("gotcha", "walltime")
    log:
        err="{OUTDIR}/logs/gotcha/{sample}_calling.err",
        out="{OUTDIR}/logs/gotcha/{sample}_calling.out"
    scripts:
        "../scripts/mutcalling.R"

rule merging: 
    input:
        mut_output=dynamic("{OUTDIR}/{sample}/gotcha/results/Split/Filtered/_rslurm_{letter}/results_0.RDS")
    ouput:
        merged="{}/{{sample}}/gotcha/results/Split/Filtered/MergedOuts/outs.collapsed.Rdata".format(OUTDIR)
    params:
        input_dir = units["gotcha"][wc.sample],
        output_dir = "{OUTDIR}/{sample}/gotcha/"
    threads: get_resource("gotcha", "threads")
    resources: 
        mem_mb=get_resource("gotcha", "mem_mb"),
        walltime=get_resource("gotcha", "walltime")
    log:
        err="{OUTDIR}/logs/gotcha/{sample}_merging.err",
        out="{OUTDIR}/logs/gotcha/{sample}_merging.out"
    scripts:
        "../scripts/mutmerging.R"
