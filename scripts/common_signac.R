library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38) #BSgenome.Mmusculus.UCSC.mm10
library(dplyr)
library(EnsDb.Hsapiens.v86) #for mm10
library(foreach)
library(GenomicRanges)
library(ggplot2)
library(harmony)
library(Signac)
library(Seurat)
library(data.table)


QCplot = function(df, metric, name.ID, color, max.y = NULL){
  df$Metric = df[, metric]
  df$NameID = df[, name.ID]
  ggplot(data = df, mapping = aes(y = Metric, x = NameID)) +
          geom_violin(size = 0.8, fill = color, color = "black") +
          stat_summary(fun.y = median, geom = "point", size = 5, fill = "snow", shape = 23) +
          ylab(metric) +
          xlab("") +
          ylim(c(0, ifelse(is.null(max.y), max(df$Metric), max.y))) +
          theme(text = element_text(colour = "black", size = 15), 
                axis.line = element_line(color = "black", size = 0.8),
                axis.text.x = element_text(colour = "black", size = 15, angle = 30, hjust = 1),
                axis.text.y = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", size = 0.8)
  )
}

FindCommonPeaks_scATACseq = function(samples.ID, dir.data.samples, dir.mito, MaxPeakWidth, MinPeakWidth, MinCounts, 
                                     dir.output,  granges.ens, genome, analysis.ID, name.ID, saveRDS = F, name.rds, mito = T){
  
  # Create a combined peak set
  ## Create a list of peak files
  # dir.data.samples = setNames(paste0(dir.data, samples.ID), samples.ID)
  peak.list = sapply(dir.data.samples, function(y) list(read.table(
    file = list.files(y, recursive = T, full.names = T)[grep(pattern = "filtered_peak_bc_matrix/peaks.bed", x = list.files(y, recursive = T, full.names = T))],
    col.names = c("chr", "start", "end")
  )))
  names(peak.list) = samples.ID
  
  ## Convert to Granges objects
  gr.list = lapply(peak.list, function(x) makeGRangesFromDataFrame(x))
  gr.list = foreach(obj = gr.list, .combine = "c") %do% obj 
  
  ## Combines and filter low quality peaks (length)
  combined.peaks = Signac::reduce(x = gr.list) #you can also use "disjoin()"
  peak.width = width(combined.peaks)
  combined.peaks = combined.peaks[peak.width < MaxPeakWidth & peak.width > MinPeakWidth]
  
  # Create fragment objects and cell filtering
  ## Load metadata
  meta.list = sapply(dir.data.samples, function(y) list(read.table(
    file = list.files(y, recursive = T, full.names = T)[grep(pattern = "singlecell.csv", x = list.files(y, recursive = T, full.names = T))],
    stringsAsFactors = F,
    sep = ",",
    header = T,
    row.names = 1
  )[-1,]))
  names(meta.list) = samples.ID
  
  ## Perform an initial filtering (low count cells)
  meta.list = lapply(meta.list, function(x) x[x[, "passed_filters"] > MinCounts,])
  
  ## Create fragment objects
  frag.list = sapply(samples.ID, function(y)  CreateFragmentObject(
    path = list.files(dir.data.samples[y], recursive = T, full.names = T)[grep(pattern = "fragments.tsv.gz", x = list.files(dir.data.samples[y], recursive = T, full.names = T))][1],
    cells = rownames(meta.list[[y]])))
  
  # Quantifying peaks in each dataset
  counts.list = sapply(samples.ID, function(y) FeatureMatrix(
    fragments = frag.list[[y]],
    features = combined.peaks,
    cells = rownames(meta.list[[y]])))
  
  # Create chromatin assay
  chrom.list = sapply(samples.ID, function(y) CreateChromatinAssay(counts.list[[y]], fragments = frag.list[[y]]))
  
  # Create Seurat objects
  seurat.list = sapply(samples.ID, function(y) CreateSeuratObject(chrom.list[[y]], assay = "ATAC", meta.data = meta.list[[y]]))
  seurat.list = sapply(samples.ID, function(y) AddMetaData(seurat.list[[y]], metadata = y, col.name = name.ID))
  
  # Load Genomic annotations
  gene_annotation = GetGRangesFromEnsDb(ensdb = granges.ens)
  gene_annotation = renameSeqlevels(gene_annotation, mapSeqlevels(seqlevels(gene_annotation), "UCSC"))
  genome(gene_annotation) = genome
  
  # Add genomic annotations
  for(y in 1:length(seurat.list)){
    Annotation(seurat.list[[y]]) = gene_annotation
  }
  
  # Compute ATAC QC metrics
  seurat.list = sapply(samples.ID, function(y) NucleosomeSignal(object = seurat.list[[y]]))
  seurat.list = sapply(samples.ID, function(y) TSSEnrichment(object = seurat.list[[y]], fast = FALSE))
  seurat.list = sapply(samples.ID, function(y) AddMetaData(seurat.list[[y]], metadata = seurat.list[[y]][["peak_region_fragments"]]/seurat.list[[y]][["passed_filters"]]*100, col.name = "pct_reads_in_peaks"))
  seurat.list = sapply(samples.ID, function(y) AddMetaData(seurat.list[[y]], metadata = seurat.list[[y]][["blacklist_region_fragments"]]/seurat.list[[y]][["peak_region_fragments"]], col.name = "blacklist_ratio"))
  
  # Add mito depth info if mito == T
  if(mito){
    mgatk.list = list()
    for (y in samples.ID){
      mgatk = ReadMGATK(dir.mito[y])
      mgatk.list[[y]] = mgatk #Will be needed further in another function
      mgatk.assay = CreateAssayObject(counts = mgatk$counts)
      mgatk.assay = subset(mgatk.assay, cells = intersect(Cells(seurat.list[[y]]), Cells(mgatk.assay)))
      seurat.list[[y]] = subset(seurat.list[[y]], cells = Cells(mgatk.assay))
      seurat.list[[y]][["mito"]] = mgatk.assay
      mgatk.depth = mgatk$depth
      seurat.list[[y]] = AddMetaData(seurat.list[[y]], mgatk.depth, col.name = "mtDNA_depth")
    }
  }
  
  # Collect all metadata
  if(saveRDS){
    dir.create(paste0(dir.output, "Objects/"), recursive = T)
  }
  
  df.meta = foreach(y = samples.ID, .combine = "rbind") %do% {
    if(saveRDS){
      saveRDS(seurat.list[[y]], paste0(dir.output, "Objects/", name.rds, "_", y, "_", analysis.ID, ".rds"))
    }
    seurat.list[[y]]@meta.data
  }
  
  # Quality control plots
  
  dir.create(paste0(dir.output, "Plots/QC/"), recursive = T)
  pdf(paste0(dir.output, "Plots/QC/QCbefore-MetricsCombined_", analysis.ID, ".pdf"), width = (4+1.3*length(samples.ID)), height = 4)
  print(QCplot(df.meta, 'pct_reads_in_peaks', name.ID, "#006600"))
  print(QCplot(df.meta, 'peak_region_fragments', name.ID, "#006666", 30000))
  print(QCplot(df.meta, 'TSS.enrichment', name.ID, "#99CC99", 50))
  print(QCplot(df.meta, 'nucleosome_signal', name.ID, "#CCCC33", 8))
  if(mito){
    print(QCplot(df.meta, 'mtDNA_depth', name.ID, "#CC9900"))
  }
  dev.off()
  
  
  return(list(Seurat = seurat.list, MGATK = mgatk.list))
}

MakeDoubletList_scATACseq = function(samples.ID, dir.doublets){
  doublets.list = lapply(1:length(samples.ID), function(y) read.table(dir.doublets[y])$V1)
  names(doublets.list) = samples.ID
  
  return(doublets.list)
}

FilterCells_scATACseq = function(seurat.obj, samples.ID, analysis.ID, dir.output, doublets, saveRDS = F, name.rds, minPeakFrag, maxPeakFrag,
                                 minPercPeak, maxNucleoSig, minTSS, minMitoDepth, mito = F){
  
  if(mito){
      seurat.obj = subset(x = seurat.obj, subset = peak_region_fragments > minPeakFrag &
               peak_region_fragments < maxPeakFrag &
               pct_reads_in_peaks > minPercPeak &
               nucleosome_signal < maxNucleoSig &
               TSS.enrichment > minTSS &
               mtDNA_depth > minMitoDepth
      )
  } else {
      seurat.obj = subset(x = seurat.obj, subset = peak_region_fragments > minPeakFrag &
               peak_region_fragments < maxPeakFrag &
               pct_reads_in_peaks > minPercPeak &
               nucleosome_signal < maxNucleoSig &
               TSS.enrichment > minTSS
      )
  }
  
  seurat.obj = subset(x = seurat.obj, cells = setdiff(Cells(seurat.obj), doublets))

  if(saveRDS){
    dir.create(paste0(dir.output, "Objects/"), recursive = T)
    saveRDS(seurat.obj, paste0(dir.output, "Objects/", name.rds, "_", samples.ID, "_", analysis.ID, ".rds"))

  }
  
  return(seurat.obj)
}

IdentVariant_scATACseq = function(seurat.obj, mgatk, samples.ID, analysis.ID, dir.output, saveVariant, MinCellVar){
    
    # Find informative mtDNA variants
    dir.create(paste0(dir.output, "Plots/Variant/"), recursive = T)
    variable.sites = IdentifyVariants(seurat.obj, assay = "mito", refallele = mgatk$refallele)
    
    pdf(paste0(dir.output, "Plots/Variant/VariantPlot_", samples.ID, "_", analysis.ID, ".pdf"), width = 6, height = 6)
    print(VariantPlot(variants = variable.sites, min.cells = MinCellVar, concordance.threshold = 1, vmr.threshold = 1) +
            theme(text = element_text(colour = "black", size = 15),
                  axis.line = element_line(color = "black", size = 0.8),
                  axis.text.x = element_text(colour = "black", size = 15, angle = 30, hjust = 1),
                  axis.text.y = element_text(colour = "black", size = 15),
                  axis.ticks = element_line(color = "black", size = 0.8)))
    dev.off()
    
    if(saveVariant){
      dir.create(paste0(dir.output, "Tables/Variants/"), recursive = T)
      fwrite(variable.sites, paste0(dir.output, "Tables/Variants/VariantTable_", samples.ID, "_", analysis.ID, ".txt"), sep = "\t")
      
    }
    
    return(variable.sites)
}

FindClonotypes_scATACseq = function(seurat.obj, variable.sites, samples.ID, dir.output, saveRDS, name.rds, analysis.ID, MinStrandCor, MinCellVar, MinVMR){
  
  # Find informative mtDNA variants
  dir.create(paste0(dir.output, "Plots/Variant/"), recursive = T)

  pdf(paste0(dir.output, "Plots/Variant/VariantPlot_Threshold_", samples.ID, "_", analysis.ID, ".pdf"), width = 6, height = 6)
  print(VariantPlot(variants = variable.sites, min.cells = MinCellVar, concordance.threshold = MinStrandCor, vmr.threshold = MinVMR) +
          theme(text = element_text(colour = "black", size = 15),
                axis.line = element_line(color = "black", size = 0.8),
                axis.text.x = element_text(colour = "black", size = 15, angle = 30, hjust = 1),
                axis.text.y = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  high.conf = subset(
    variable.sites, subset = n_cells_conf_detected >= MinCellVar &
      strand_correlation >= MinStrandCor &
      vmr > MinVMR
  )
  
  # Compute the variant allele frequency for each cell
  seurat.obj = AlleleFreq(
    object = seurat.obj,
    variants = high.conf$variant,
    assay = "mito"
  )
  
  DefaultAssay(seurat.obj) = "alleles"
  seurat.obj = FindClonotypes(seurat.obj)
  
  if(saveRDS){
    dir.create(paste0(dir.output, "Objects/"), recursive = T)
    saveRDS(seurat.obj, paste0(dir.output, "Objects/", name.rds, "_", samples.ID, "_", analysis.ID, ".rds"))
    
  }
  
  return(seurat.obj)
}

Merge_scATACseq = function(seurat.list, dir.output, saveRDS = F, name.rds, analysis.ID, name.ID, mito = T){
  samples.ID = names(seurat.list)
  
  combined = merge(
    x = seurat.list[[1]],
    y = seurat.list[names(seurat.list) != samples.ID[1]],
    add.cell.ids = samples.ID
  )
  
  pdf(paste0(dir.output, "Plots/QC/QCafter-Metrics_", analysis.ID, ".pdf"), width = (4+1.3*length(samples.ID)), height = 4)
  print(QCplot(combined@meta.data, 'pct_reads_in_peaks', name.ID, "#006600"))
  print(QCplot(combined@meta.data, 'peak_region_fragments', name.ID, "#006666", 30000))
  print(QCplot(combined@meta.data, 'TSS.enrichment', name.ID, "#99CC99", 50))
  print(QCplot(combined@meta.data, 'nucleosome_signal', name.ID, "#CCCC33", 8))
  if(mito){
    print(QCplot(combined@meta.data, 'mtDNA_depth', name.ID, "#CC9900"))
  }
  dev.off()
  
  if(saveRDS){
    dir.create(paste0(dir.output, "Objects/"), recursive = T)
    saveRDS(combined, paste0(dir.output, "Objects/", name.rds, "_", analysis.ID, ".rds"))
  }
  return(combined)
}

GeneActivity_scATACseq = function(seurat.obj, saveRDS = F, name.rds, analysis.ID, dir.output){
  DefaultAssay(seurat.obj) = "ATAC"
  
  # Compute gene accessibility
  gene.activities = GeneActivity(seurat.obj)
  
  # Add to the Seurat object as a new assay and normalize
  seurat.obj[['RNA']] = CreateAssayObject(counts = gene.activities)
  
  seurat.obj = NormalizeData(
    object = seurat.obj,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(seurat.obj$nCount_RNA)
  )
  
  if(saveRDS){
    dir.create(paste0(dir.output, "Objects/"), recursive = T)
    saveRDS(seurat.obj, paste0(dir.output, "Objects/", name.rds, "_", analysis.ID, ".rds"))
  }
  
  return(seurat.obj)
}
