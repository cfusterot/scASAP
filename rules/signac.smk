import glob

if config["signac"]["enabled"]:
    rule signac:
        input:
            h5="{}/{{sample}}/cellranger_count/outs/filtered_peak_bc_matrix.h5".format(OUTDIR),
            fragpath="{}/{{sample}}/cellranger_count/outs/fragments.tsv.gz".format(OUTDIR), 
            metadata="{}/{{sample}}/cellranger_count/outs/singlecell.csv".format(OUTDIR),
            mgatk="{}/{{sample}}/mgatk/final/{{sample}}.variant_stats.tsv.gz".format(OUTDIR),
            amulet="{}/{{sample}}/amulet/MultipletSummary.txt".format(OUTDIR)
        output:
            integrated="{}/{{sample}}/signac/integrated_preqc.rds".format(OUTDIR)
        conda:
            "../envs/signac.yaml"
        params:
            min_peak_width=config['signac']['min_peak_width'],
            max_peak_width=config['signac']['max_peak_width'],
            min_counts=config['signac']['min_counts'],
            max_nucleosome_signal=config['signac']['max_nucleosome_signal'],
            min_TSS=config['signac']['min_TSS'],
            min_peak_fragment=config['signac']['min_peak_fragment'],
            max_peak_fragment=config['signac']['max_peak_fragment'],
            min_percentage_peaks = config['signac']["min_percentage_peaks"],
            min_depth=config['signac']['mito']['min_depth'],
            min_strand_cor=config['signac']['mgatk']['min_strand_cor'],
            min_cell_var=config['signac']['mgatk']['min_cell_var'],
            min_VMR=config['signac']['mgatk']['min_VMR']
        threads: get_resource("signac", "threads")
        resources:
            mem_mb=get_resource("signac", "mem_mb"),
            walltime=get_resource("signac", "walltime")
        log:
            err="{}/{{sample}}/signac.err".format(LOGDIR),
            out="{}/{{sample}}/signac.out".format(LOGDIR)
        # Add here something 
        script:
            "../scripts/signac.R"
