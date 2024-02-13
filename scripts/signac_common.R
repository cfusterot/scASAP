# ------ Signac Analysis functions ------ #


# --- 1. Find common set of peaks in a collection of samples --- #
# Parameter list:
# - dir.data.samples
# - MaxPeakWidth
# - MinPeakWidth
# - dir.output 
FindCommonPeaks_scATACseq = function(dir.data.samples,  MaxPeakWidth, MinPeakWidth, dir.output){
  # Create a combined peak set
  ## Create a list of peak files
  peak.list = sapply(dir.data.samples, function(y)
                     list(read.table(paste0(y,"/filtered_peak_bc_matrix/peaks.bed"),
                                     col.names = c("chr", "start", "end"))))
  ## Convert to Granges objects
  gr.list = lapply(peak.list, function(x) makeGRangesFromDataFrame(x))
  gr.list = foreach(obj = gr.list, .combine = "c") %do% obj 
                       
  ## Combines and filter low quality peaks (length)
  combined.peaks = Signac::reduce(x = gr.list) #you can also use "disjoin()"
  peak.width = width(combined.peaks)
  combined.peaks = combined.peaks[peak.width < MaxPeakWidth & peak.width > MinPeakWidth]
  write.table(combined.peaks, paste0(dir.output, "/CommonSetOfPeaks.bed"))
 }

# --- 2. Generate ATAC object  --- #
# Parameter list:
# - sample.ID
# - dir.data.sample
# - MinCounts
# - dir.output
# - gannotation
# - name.ID
CreateATACobject = function(sample.ID, dir.data.sample, MinCounts, dir.output, gannotation, name.ID){
    if(integration){
          combined.peaks = read.table(paste0(dir.output, "/CommonSetOfPeaks.bed"))
  } else {
        combined.peaks = read.table(paste0(dir.data.sample, "/filtered_peak_bc_matrix/peaks.bed"))
    }
  combined.peaks = makeGRangesFromDataFrame(combined.peaks)
  # Create fragment objects and cell filtering
  # Load metadata
  meta = read.table(
    file = paste0(dir.data.sample, "/cellranger_count/singlecell.csv"),
    stringsAsFactors = F,
    sep = ",",
    header = T,
    row.names = 1
  )[-1,]
  # Perform an initial filtering (low count cells)
  meta = meta[meta[, "passed_filters"] > MinCounts,]
  cells.to.keep = read.table(paste0(dir.data.sample, "/cellranger_count/filtered_peak_bc_matrix/barcodes.tsv"))$V1 # Discard empty droplets
  meta = meta[cells.to.keep,]                        
  # Create fragment objects
  frag = CreateFragmentObject(path = paste0(dir.data.sample, "/cellranger_count/fragments.tsv.gz"), cells = rownames(meta))
  # Quantifying peaks
  counts = FeatureMatrix(fragments = frag, features = combined.peaks, cells = rownames(meta))
  # Create chromatin assay
  chrom = CreateChromatinAssay(counts, fragments = frag)              
  # Create Seurat objects
  chrom = CreateSeuratObject(chrom, assay = "peaks", meta.data = meta)
  chrom = AddMetaData(chrom, metadata = sample.ID, col.name = name.ID)
  # Remove doublets
  doublets = read.table(paste0(dir.data.sample, "/amulet/MultipletBarcodes_01.txt"))$V1
  chrom = subset(chrom, cells = setdiff(Cells(chrom), doublets))
  # Add genomic annotations
  Annotation(chrom) = gannotation                                                                 
  # Compute ATAC QC metrics
  chrom = NucleosomeSignal(object = chrom)
  chrom = TSSEnrichment(object = chrom, fast = FALSE)
  chrom = AddMetaData(chrom, metadata = chrom[["peak_region_fragments"]]/chrom[["passed_filters"]]*100, 
                      col.name = "pct_reads_in_peaks")
  saveRDS(chrom, paste0(dir.output, "/SeuratObject_", sample.ID,".rds"))
  # Collect all metadata
  return(chrom@meta.data)
}

# --- 3. QC plot generation  --- #
# Parameter list:
# - df
# - metric 
# - color
# - n.samples
# - name.ID 
# - max.y = NULL
QCplot = function(df, metric, color, n.samples, name.ID, max.y = NULL){
  df$Metric = df[, metric]
  df$NameID = df[, name.ID]
  ggplot(data = df, mapping = aes(y = Metric, x = NameID)) +
    geom_violin(linewidth = 0.8, fill = color, color = "black") +
    stat_summary(fun.y = median, geom = "point", size = 3, fill = "snow", shape = 23) +
    ylab(metric) + xlab("") + ylim(c(0, ifelse(is.null(max.y), max(df$Metric), max.y))) +
    theme_classic() + theme(text = element_text(colour = "black", size = 15), 
                            aspect.ratio = 1/(0.3*n.samples),
                            axis.line = element_line(color = "black", size = 0.8),
                            axis.text.x = element_text(colour = "black", size = 15, angle = 30, hjust = 1),                                                
                            axis.text.y = element_text(colour = "black", size = 15), 
                            axis.ticks = element_line(color = "black", size = 0.8))
}


# --- 4. Filter ATAC object  --- #
# Parameter list:
# - seurat
# - minPeakFrag
# - maxPeakFrag
# - minPercPeak
# - maxNucleoSig
# - minTSS
Filter_ATAC = function(seurat, minPeakFrag, maxPeakFrag, minPercPeak, maxNucleoSig, minTSS){
  seurat = subset(seurat, subset = pct_reads_in_peaks >= minPercPeak & peak_region_fragments >= minPeakFrag & peak_region_fragments < maxPeakFrag &
                   TSS.enrichment >= minTSS & nucleosome_signal < maxNucleoSig)
  return(seurat)
}

# --- 5. Add_MGATK --- #
# Parameter list:
# - seurat
# - dir.data.sample
# - sample.ID
# - MinCellVar
# - MinStrandCor 
# - MinVMR
# - minMitoDepth
Add_MGATK = function(seurat, dir.data.sample, sample.ID, MinCellVar, MinStrandCor, MinVMR, minMitoDepth){

  # Add mito depth info
  mgatk = ReadMGATK(paste0(dir.data.sample, "/mgatk/final/"))
  mgatk.assay = CreateAssayObject(counts = mgatk$counts)
  mgatk.assay = subset(mgatk.assay, cells = intersect(Cells(seurat), Cells(mgatk.assay)))
  seurat = subset(seurat, cells = Cells(mgatk.assay))
  seurat[["mito"]] = mgatk.assay
  seurat = AddMetaData(seurat, mgatk$depth, col.name = "mtDNA_depth")

  seurat = subset(seurat, subset = mtDNA_depth >= minMitoDepth)

  #Identify variants
  variable.sites = IdentifyVariants(seurat, assay = "mito", refallele = mgatk$refallele)
  dir.create(file.path(dir.data.sample, "signac", "plots"))
  pdf(paste0(dir.data.sample, "/signac/plots/VariantPlot_", sample.ID, ".pdf"), width = 6, height = 6)
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

  write.table(high.conf, paste0(dir.data.sample, "/signac/HighConfVariants_", sample.ID, ".txt"))
  return(seurat)
}
# --- 5. Heteroplasmy level calculation --- #
# Parameters list:
# - sample.ID 
# - variants 
# - dir.output
Calculate_Hetroplasmy = function(sample.ID, variants, dir.output){
  seurat = readRDS(paste0(dir.output, "/SeuratObject_", sample.ID, ".s2.rds"))

  seurat = AlleleFreq(
    object = seurat,
    variants = variants,
    assay = "mito"
  )

  saveRDS(seurat, paste0(dir.output, "/SeuratObject_", sample.ID, ".s2.rds"))
}

##from Caleb
seuratSNN_cosineDistance <- function(obj, resolution){
  clusters <- FindClusters(object = obj$snn, resolution = resolution)
  return(as.character(clusters[,1]))
}

getNN <- function(mat_af, k.param = 10){
  set.seed(1)
  rownames(mat_af) <- make.unique(rownames(mat_af))
  obj <- FindNeighbors(mat_af, k.param = k.param, annoy.metric = "cosine")
  obj
}

# --- 6. Final Signac workflow --- #
# Parameter list
# - seurat
# - dir.output
# - name.grp
# - harmony
# - CutOff_FTF
# - NbDim
# - RmCompo
# - Reduction
# - MinDistUMAP
# - AlgoClustering
# - ResolutionClustering
# - col.palette
SignacWorkflow = function(seurat, dir.output, name.grp, harmony, CutOff_FTF, NbDim, RmCompo, Reduction,  MinDistUMAP, AlgoClustering, ResolutionClustering, col.palette) {
  DefaultAssay(seurat) = "peaks"

  if(!harmony){
    seurat = RunTFIDF(seurat)
    seurat = FindTopFeatures(seurat, min.cutoff = CutOff_FTF)
    seurat = RunSVD(seurat)
    
    pdf(paste0(dir.output, "/Plots/DepthCor.pdf"), width = 6, height = 4)
    print(DepthCor(seurat) + theme(text = element_text(colour = "black", size = 15), 
                                   axis.line = element_line(color = "black", size = 0.8),
                                   axis.text.x = element_text(colour = "black", size = 15),
                                   axis.text.y = element_text(colour = "black", size = 15),
                                   axis.ticks = element_line(color = "black", size = 0.8)))
    dev.off()
  }

  
  if(RmCompo == TRUE){
    Dim = 2:NbDim
  } else if(RmCompo == FALSE){
    Dim = 1:NbDim
  }
  
  seurat = RunUMAP(seurat, dims = Dim, reduction = Reduction, min.dist = MinDistUMAP)
  
  pdf(paste0(dir.output, "/Plots/UMAP_Merge_Condition_", Reduction, ".pdf"), width = 6, height = 4)
  sapply(name.grp, function(x){
    print(DimPlot(seurat, group.by = x, pt.size = 0.001, shuffle = T, cols = col.palette) + 
            theme(text = element_text(colour = "black", size = 15), 
                  axis.line = element_line(color = "black", size = 0.8),
                  axis.text.x = element_text(colour = "black", size = 15),
                  axis.text.y = element_text(colour = "black", size = 15),
                  axis.ticks = element_line(color = "black", size = 0.8)))
  })
  dev.off()
  
  seurat = FindNeighbors(object = seurat, reduction = Reduction, dims = Dim)
  seurat = FindClusters(object = seurat, verbose = FALSE, algorithm = AlgoClustering, resolution = ResolutionClustering)
  
  pdf(paste0(dir.output, "/Plots/UMAP_Merge_Clusters_", Reduction, ".pdf"), width = 6, height = 4)
  print(DimPlot(seurat, group.by = "seurat_clusters", pt.size = 0.001, shuffle = T, cols = col.palette) + 
            theme(text = element_text(colour = "black", size = 15), 
                  axis.line = element_line(color = "black", size = 0.8),
                  axis.text.x = element_text(colour = "black", size = 15),
                  axis.text.y = element_text(colour = "black", size = 15),
                  axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()

  return(seurat)
}
