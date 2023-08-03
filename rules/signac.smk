import glob

if config["signac"]["enable"]:
    rule step1_signac:
        input:
            outs=get_integration_outs, 
            mgatk=get_integration_mgatk, 
            amulet=get_integration_amulet
        output:
            directory=expand("{OUTDIR}/signac", OUTDIR = OUTDIR),
            qc_plots=expand("{OUTDIR}/signac/plots/vlnplot_qc_beforefiltering_ATAC.pdf", OUTDIR = OUTDIR)
        conda:
            "../envs/signac.yaml"
        params: 
            min_peak_width=config['signac']['qc']['min_peak_width'],
            max_peak_width=config['signac']['qc']['max_peak_width'],
            min_counts=config['signac']['qc']['min_counts'],
            integration=config['signac']['integrate_samples'],
            nb_cores=config['signac']['cores'],
            granges_ens=config['signac']['annotation']['granges_ens'],
            genome=config['signac']['annotation']['genome']
        threads: get_resource("signac", "threads")
        resources:
            mem_mb=get_resource("signac", "mem_mb"),
            walltime=get_resource("signac", "walltime")
        log:
            err=expand("{LOGDIR}/signac/step1_signac.err", LOGDIR = LOGDIR),
            out=expand("{LOGDIR}/signac/step1_signac.out", LOGDIR = LOGDIR)
        # Add here something 
        script:
            "../scripts/step1_qc.R"

    rule step2_signac:
        input:
            seurat_in = expand("{OUTDIR}/signac/SeuratObject_{alias}.rds", OUTDIR=OUTDIR, alias=alias)
        output:
            seurat_out = expand("{OUTDIR}/signac/SeuratObject_{alias}.rds", OUTDIR=OUTDIR, alias=alias) 
        conda:
            "../envs/signac.yaml"
        params:
            mito = config['signac']['mito']['enable'],
            integration = config['signac']['integrate_samples'],
            cores = config['signac']['cores'],
            min_peak_fragment = config['signac']['qc']['min_peak_fragment'],
            max_peak_fragment = config['signac']['qc']['max_peak_fragment'],
            min_percentage_peaks = config['signac']['qc']['min_percentage_peaks'], 
            max_nucleosome_signal = config['signac']['qc']['max_nucleosome_signal'],
            min_TSS = config['signac']['qc']['min_TSS'],
            min_strand_cor = config['signac']['qc']['min_strand_cor'],
            min_VMR = config['signac']['mito']['min_VMR'],
            min_cell_var = config['signac']['mito']['min_cell_var'],
            min_depth = config['signac']['mito']['min_depth']
        threads: get_resource("signac", "threads")
        resources:
            mem_mb=get_resource("signac", "mem_mb"),
            walltime=get_resource("signac", "walltime")
        log:
            err=expand("{LOGDIR}/signac/step2_signac.err", LOGDIR = LOGDIR),
            out=expand("{LOGDIR}/signac/step2_signac.out", LOGDIR = LOGDIR)
        script:
            "../scripts/step2_mitoassay.R"

    rule step3_signac:
        input:
            seurat_in = expand("{OUTDIR}/signac/SeuratObject_{alias}.rds", OUTDIR=OUTDIR, alias=alias)
        output:
            seurat_out = expand("{OUTDIR}/signac/SeuratObject_{condition}.rds", OUTDIR=OUTDIR, condition=condition)
        conda:
            "../envs/signac.yaml"
        params:
            cores = config['signac']['cores'],
            fc_resolution = config['signac']['mito']['clonotype_finding']['fc_resolution'],
            fc_k = config['signac']['mito']['clonotype_finding']['fc_k']
        threads: get_resource("signac", "threads")
        resources:
        log:
            err=expand("{LOGDIR}/signac/step2_signac.err", LOGDIR = LOGDIR),
            out=expand("{LOGDIR}/signac/step2_signac.out", LOGDIR = LOGDIR)
        script:
            "../scripts/step3_findclonotypes.R"

    rule step4_signac:
        input:
            seurat_out = expand("{OUTDIR}/signac/SeuratObject_{condition}.rds", OUTDIR=OUTDIR, condition=condition)
        output:
            expand("{OUTDIR}/signac/SeuratObject_Merge.rds", OUTDIR = OUTDIR)
        conda:
            "../envs/signac.yaml"
        params:
           cores = config['signac']['cores'],
           harmony = config['signac']['harmony'],
           lambdaHarmony = config['signac']['harmony']['lambda'],
           GEX = config['signac']['GEX'],
           RmCompo = config['signac']['clustering']['remove_first_component'],
           NbDim = config['signac']['clustering']['dims'],
           CutOff_FTF = config['signac']['clustering']['min.cutoff'],
           MinDistUMAP = config['signac']['clustering']['min.dist'],
           AlgoClustering = config['signac']['clustering']['algorithm'],
           ResolutionClustering = config['signac']['clustering']['resolution'],
           batch_corr = config['signac']['clustering']['batch_correction'],
           var.batch  = config['signac']['clustering']['batch_correction']['variable_to_regress'],
           nb.cores = config['signac']['cores']
        threads: get_resource("signac", "threads")
        resources:
        log:
            err=expand("{LOGDIR}/signac/step2_signac.err", LOGDIR = LOGDIR),
            out=expand("{LOGDIR}/signac/step2_signac.out", LOGDIR = LOGDIR)
        script:
            "../scripts/step4_merge.R"
