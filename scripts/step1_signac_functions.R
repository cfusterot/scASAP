FindCommonPeaks_scATACseq = function(dir.data.samples,  MaxPeakWidth, MinPeakWidth, dir.output){
  # Create a combined peak set
  ## Create a list of peak files
  peak.list = sapply(dir.data.samples, function(y)
                     list(read.table(paste0(y,"/cellranger_count/outs/filtered_peak_bc_matrix/peaks.bed"),
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

CreateATACobject = function(sample.ID, dir.data.sample, MinCounts, dir.output, gannotation, name.ID){
    if(integration){
          combined.peaks = read.table(paste0(dir.output, "/CommonSetOfPeaks.bed"))
  } else {
        combined.peaks = read.table(paste0(dir.data.sample, "/cellranger_count/outs/filtered_peak_bc_matrix/peaks.bed"))
    }
  combined.peaks = makeGRangesFromDataFrame(combined.peaks)
  # Create fragment objects and cell filtering
  # Load metadata
  meta = read.table(
    file = paste0(dir.data.sample, "/cellranger_count/outs/singlecell.csv"),
    stringsAsFactors = F,
    sep = ",",
    header = T,
    row.names = 1
  )[-1,]
  # Perform an initial filtering (low count cells)
  meta = meta[meta[, "passed_filters"] > MinCounts,]
  cells.to.keep = read.table(paste0(dir.data.sample, "/cellranger_count/outs/filtered_peak_bc_matrix/barcodes.tsv"))$V1 # Discard empty droplets
  meta = meta[cells.to.keep,]                        
  # Create fragment objects
  frag = CreateFragmentObject(path = paste0(dir.data.sample, "/cellranger_count/outs/fragments.tsv.gz"), cells = rownames(meta))
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
                      col.name = "pct_reads_in_peaks")                                                                           saveRDS(chrom, paste0(dir.output, "/SeuratObject_", sample.ID,".rds"))
  # Collect all metadata
  return(chrom@meta.data)
}

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
                            axis.text.x = element_text(colour = "black", size = 15, angle = 30, hjust = 1),                                                axis.text.y = element_text(colour = "black", size = 15), 
                            axis.ticks = element_line(color = "black", size = 0.8))
}
