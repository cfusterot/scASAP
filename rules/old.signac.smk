import glob

if config["signac"]["enable"]:
    rule signac:
        input:
            outs = get_integration_outs, 
            mgatk = get_integration_mgatk, 
            amulet = get_integration_amulet
        output:
            directory=expand("{OUTDIR}/signac", OUTDIR = OUTDIR),
            integrated=expand("{OUTDIR}/signac/SeuratCombined_GeneScore.rds", OUTDIR = OUTDIR)
        conda:
            "../envs/signac.yaml"
        params:
            min_peak_width=config['signac']['qc']['min_peak_width'],
            max_peak_width=config['signac']['qc']['max_peak_width'],
            min_counts=config['signac']['qc']['min_counts'],
            max_nucleosome_signal=config['signac']['qc']['max_nucleosome_signal'],
            min_TSS=config['signac']['qc']['min_TSS'],
            min_peak_fragment=config['signac']['qc']['min_peak_fragment'],
            max_peak_fragment=config['signac']['qc']['max_peak_fragment'],
            min_percentage_peaks = config['signac']['qc']["min_percentage_peaks"],
            min_depth=config['signac']['mito']['min_depth'],
            min_strand_cor=config['signac']['mito']['min_strand_cor'],
            min_cell_var=config['signac']['mito']['min_cell_var'],
            min_VMR=config['signac']['mito']['min_VMR']
        threads: get_resource("signac", "threads")
        resources:
            mem_mb=get_resource("signac", "mem_mb"),
            walltime=get_resource("signac", "walltime")
        log:
            err=expand("{LOGDIR}/signac/signac.err", LOGDIR = LOGDIR),
            out=expand("{LOGDIR}/signac/signac.out", LOGDIR = LOGDIR)
        # Add here something 
        script:
            "../scripts/signac.R"
