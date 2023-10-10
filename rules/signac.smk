import glob

if config["signac"]["enable"]:
    if config["signac"]["integrate_samples"]:
        rule step0_peaks:
            input:
                get_common_peaks
            output:
                peaks="{}/integration/CommonSetOfPeaks.bed".format(OUTDIR)
            conda:
                "../envs/signac.yaml"
            params:
                directory=directory("{}/integration/").format(OUTDIR),
                min_peak_width=config['signac']['qc']['min_peak_width'],
                max_peak_width=config['signac']['qc']['max_peak_width']
            threads:  get_resource("signac", "threads")
            resources:
                mem_mb=get_resource("signac", "mem_mb"),
                walltime=get_resource("signac", "walltime")
            log:
                err="{}/signac/step0_peaks.err".format(LOGDIR),
                out="{}/signac/step0_peaks.out".format(LOGDIR)
            script:
                "../scripts/step0_peaks.R"
    
        rule step1_signac:
            input:
                peaks="{}/integration/CommonSetOfPeaks.bed".format(OUTDIR),
                out="{}/{{sample}}/cellranger_count/cellranger.finish".format(OUTDIR),
                mgatk="{}/{{sample}}/mgatk/final/{{sample}}.variant_stats.tsv.gz".format(OUTDIR),
                amulet="{}/{{sample}}/amulet/MultipletSummary.txt".format(OUTDIR)
            output:
                file="{}/{{sample}}/signac/SeuratObject_{{sample}}.rds".format(OUTDIR)
            conda:
                "../envs/signac.yaml"
            params: 
                sample_ID="{{sample}}",
                directory=directory("{}/{{sample}}/signac/").format(OUTDIR),
                mgatk_dir=directory("{}/{{sample}}/mgatk/final/").format(OUTDIR),
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
                err="{}/signac/{{sample}}/step1_signac.err".format(LOGDIR),
                out="{}/signac/{{sample}}/step1_signac.out".format(LOGDIR)
            script:
                "../scripts/step1_qc.R"

        rule step2_signac:
            input:
                file = "{}/{{alias}}/signac/SeuratObject_{{alias}}.rds".format(OUTDIR)
            output:
                file = "{}/{{alias}}/signac/SeuratObject_{{alias}}.s2.rds".format(OUTDIR)
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
                min_strand_cor = config['signac']['mito']['min_strand_cor'],
                min_VMR = config['signac']['mito']['min_VMR'],
                min_cell_var = config['signac']['mito']['min_cell_var'],
                min_depth = config['signac']['mito']['min_depth']
            threads: get_resource("signac", "threads")
            resources:
                mem_mb=get_resource("signac", "mem_mb"),
                walltime=get_resource("signac", "walltime")
            log:
                err="{}/signac/{{alias}}_step2_signac.err".format(LOGDIR),
                out="{}/signac/{{alias}}_step2_signac.out".format(LOGDIR)
            script:
                "../scripts/step2_mitoassay.R"

        rule step3_signac:
            input:
                seurat_in = "{}/{{alias}}/signac/SeuratObject_{{alias}}.s2.rds".format(OUTDIR)
            output:
                seurat_out = "{}/{{alias}}/signac/SeuratObjectBis_{{alias}}.rds".format(OUTDIR)
            conda:
                "../envs/signac.yaml"
            params:
                fc_resolution = config['signac']['mito']['clonotype_finding']['fc_resolution'],
                fc_k = config['signac']['mito']['clonotype_finding']['fc_k']
            threads: get_resource("signac", "threads")
            resources:
            log:
                err="{}/signac/{{alias}}_step2_signac.err".format(LOGDIR),
                out="{}/signac/{{alias}}_step2_signac.out".format(LOGDIR)
            script:
                "../scripts/step3_findclonotypes.R"

        rule step4_signac:
            input:
                get_step3_output
            output:
                "{}/integration/SeuratObject_Merge.rds".format(OUTDIR)
            conda:
                "../envs/signac.yaml"
            params:
                cores = config['signac']['cores'],
                harmony = config['signac']['harmony']['enable'],
                lambdaHarmony = config['signac']['harmony']['lambda'],
                GEX = config['signac']['GEX'],
                RmCompo = config['signac']['clustering']['remove_first_component'],
                NbDim = config['signac']['clustering']['dims'],
                CutOff_FTF = config['signac']['clustering']['min.cutoff'],
                MinDistUMAP = config['signac']['clustering']['min.dist'],
                AlgoClustering = config['signac']['clustering']['algorithm'],
                ResolutionClustering = config['signac']['clustering']['resolution'],
                batch_corr = config['signac']['batch_correction']['enable'],
                var_batch  = config['signac']['batch_correction']['variable_to_regress'],
                nb_cores = config['signac']['cores']
            threads: get_resource("signac", "threads")
            resources:
            log:
                err=expand("{LOGDIR}/signac/step2_signac.err", LOGDIR = LOGDIR),
                out=expand("{LOGDIR}/signac/step2_signac.out", LOGDIR = LOGDIR)
            script:
                "../scripts/step4_merge.R"

### Triggers invidiual sample analysis ### 
    if config["signac"]["individual_analysis"]:
        rule step1_preprocess:
            input:
                outs="{}/{{sample}}/cellranger_count/outs/".format(OUTDIR), 
                mgatk="{}/{{sample}}/mgatk/final".format(OUTDIR), 
                amulet="{}/{{sample}}/amulet/MultipletBarcodes_01.txt".format(OUTDIR)
            output:
                report="{}/{{sample}}/signac/01_preprocessing_{{sample}}.html".format(OUTDIR)
            conda:
                "../envs/signac.yaml"
            params: 
                directory = expand("{OUTDIR}/", OUTDIR = OUTDIR),
                reference = config['signac']['annotation']['reference'],
                haplogroup = config['signac']['annotation']['haplogroup'],
                ncount_atac_max = config['signac']['qc']['max_counts'],
                ncount_atac_min = config['signac']['qc']['min_counts'],
                nuclesome_signal = config['signac']['qc']['max_nucleosome_signal'],
                tss_enrichment = config['signac']['qc']['min_TSS'],
                min_depth = config['signac']['mito']['min_depth'],
                n_cells_conf_detected = config['signac']['mito']['n_cells_conf_detected'],
                strand_correlation = config['signac']['mito']['min_strand_cor'],
                min_cell_var = config['signac']['mito']['min_cell_var']
            threads: get_resource("signac", "threads")
            resources:
                mem_mb=get_resource("signac", "mem_mb"),
                walltime=get_resource("signac", "walltime")
            log:
                err=expand("{LOGDIR}/signac/{{sample}}_preprocess.err", LOGDIR = LOGDIR),
                out=expand("{LOGDIR}/signac/{{sample}}_preprocess.out", LOGDIR = LOGDIR)
            script:
                "../scripts/01_preprocess_render.R"


