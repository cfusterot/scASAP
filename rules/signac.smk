import glob

if config["signac"]["enable"]:
    if config["signac"]["integrate_samples"]:
        rule step0_peaks:
            input:
                get_cellranger_finish
            output:
                peaks="{}/integration/CommonSetOfPeaks.bed".format(OUTDIR)
            conda:
                "../envs/signac.yaml"
            params:
                input_files=get_common_peaks,
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
                mgatk="{}/{{sample}}/mgatk/final/{{sample}}.rds".format(OUTDIR),
                amulet="{}/{{sample}}/amulet/MultipletSummary.txt".format(OUTDIR)
            output:
                file="{}/integration/SeuratObject_{{sample}}.rds".format(OUTDIR)
            conda:
                "../envs/signac.yaml"
            params: 
                dir_sample=directory("{}/{{sample}}/").format(OUTDIR),
                dir_integration=directory("{}/integration/").format(OUTDIR),
                sample_ID="{{sample}}",
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
                file = "{}/integration/SeuratObject_{{sample}}.rds".format(OUTDIR)
            output:
                file = "{}/{{sample}}/signac/SeuratObject_{{sample}}.s2.rds".format(OUTDIR)
            conda:
                "../envs/signac.yaml"
            params:
                dir_sample=directory("{}/{{sample}}/").format(OUTDIR),
                sample_ID="{{sample}}",
                mito = config['signac']['mito']['enable'],
                integration = config['signac']['integrate_samples'],
                dir_integration=directory("{}/integration/").format(OUTDIR),
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
                err="{}/signac/{{sample}}_step2_signac.err".format(LOGDIR),
                out="{}/signac/{{sample}}_step2_signac.out".format(LOGDIR)
            script:
                "../scripts/step2_mitoassay.R"

        rule step3_signac:
            input:
               seurat_in = lambda wildcards:expand("{outdir}/{sample}/signac/SeuratObject_{sample}.s2.rds", sample=condition_to_samples[wildcards.condition], outdir=OUTDIR)
            output:
                seurat_out = "{}/integration/SeuratObjectBis_{{condition}}.rds".format(OUTDIR)
            conda:
                "../envs/signac.yaml"
            params:
                dir_sample = lambda wildcards: expand("{outdir}/{{sample}}/", sample=condition_to_samples[wildcards.condition], outdir=OUTDIR),
                sample_ID = lambda wildcards: condition_to_samples[wildcards.condition],
                samples_tsv = config['samples'],
                fc_resolution = config['signac']['mito']['clonotype_finding']['fc_resolution'],
                fc_k = config['signac']['mito']['clonotype_finding']['fc_k']
            threads: get_resource("signac", "threads")
            resources:
                mem_mb=get_resource("signac", "mem_mb"),
                walltime=get_resource("signac", "walltime")
            log:
                err="{}/integration/{{condition}}_step3_signac.err".format(LOGDIR),
                out="{}/integration/{{condition}}_step3_signac.out".format(LOGDIR)
            script:
                "../scripts/step3_findclonotypes.R"

        #rule step4_signac:
        #    input:
        #        get_step4_input
        #    output:
        #        "{}/integration/SeuratObject_Merge.rds".format(OUTDIR)
        #    conda:
        #        "../envs/signac.yaml"
        #    params:
        #        cores = config['signac']['cores'],
        #        harmony = config['signac']['harmony']['enable'],
        #        lambdaHarmony = config['signac']['harmony']['lambda'],
        #        GEX = config['signac']['GEX'],
        #        RmCompo = config['signac']['clustering']['remove_first_component'],
        #        NbDim = config['signac']['clustering']['dims'],
        #        CutOff_FTF = config['signac']['clustering']['min.cutoff'],
        #        MinDistUMAP = config['signac']['clustering']['min.dist'],
        #        AlgoClustering = config['signac']['clustering']['algorithm'],
        #        ResolutionClustering = config['signac']['clustering']['resolution'],
        #        batch_corr = config['signac']['batch_correction']['enable'],
        #        var_batch  = config['signac']['batch_correction']['variable_to_regress'],
        #        nb_cores = config['signac']['cores']
        #    threads: get_resource("signac", "threads")
        #    resources:
        #    log:
        #        err=expand("{LOGDIR}/signac/step4_signac.err", LOGDIR = LOGDIR),
        #        out=expand("{LOGDIR}/signac/step4_signac.out", LOGDIR = LOGDIR)
        #    script:
        #        "../scripts/step4_merge.R"

### Triggers invidiual sample analysis ### 
    if config["signac"]["individual_analysis"]:
        rule step1_preprocess:
            input:
                outs="{}/{{sample}}/cellranger_count/cellranger.finish".format(OUTDIR), 
                mgatk="{}/{{sample}}/mgatk/final/{{sample}}.rds".format(OUTDIR), 
                amulet="{}/{{sample}}/amulet/MultipletSummary.txt".format(OUTDIR)
            output:
                report="{}/{{sample}}/signac/01_preprocessing_{{sample}}.html".format(OUTDIR)
            conda:
                "../envs/signac.yaml"
            params: 
                mgatk="{}/{{sample}}/mgatk/final/".format(OUTDIR),
                outs="{}/{{sample}}/cellranger_count/".format(OUTDIR),
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


