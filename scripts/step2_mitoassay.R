#################################################################
#                                                               #
#                           STEP 2:                             #
#               Filtering, add MGATK, call variants             #
#                                                               #
#################################################################
#################################################################
#                                                               #
#                  20.07.2023 - Ambre Giguelay                  #
#                                                               #
#################################################################

# -------- Load functions -------- #
message("Loading analysis functions")
source("scripts/signac_common.R")

# -------- Load libraries -------- #
library("foreach")
library("doParallel")
library("ggplot2")
library("Signac")
library("Seurat")

# -------- Load parameters -------- #
mito = snakemake@params[['mito']]
integration = snakemake@params[["integration"]]
dir.output = "processing"
name.ID = "sample"
cores = snakemake@params[["cores"]]
## chosen by user
minPeakFrag = snakemake@params[["min_peak_fragment"]]
maxPeakFrag = snakemake@params[["max_peak_fragment"]]
minPercPeak = snakemake@params[["min_percentage_peaks"]]
maxNucleoSig = snakemake@params[["max_nucleosome_signal"]]
minTSS =  snakemake@params[["min_TSS"]]
## default
MinStrandCor = snakemake@params[["min_strand_cor"]]
MinVMR = snakemake@params[["min_VMR"]]
MinCellVar = snakemake@params[["min_cell_var"]]
minMitoDepth = snakemake@params[["min_depth"]]

# -------- Read files, set variables -------- #
samples = read.table("config/samples.tsv", header = T)

# -------- Run functions -------- #
if(cores > 1){
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
} else {
  registerDoSEQ()
}

meta = foreach(j = 1:nrow(samples), .packages = c("Signac", "Seurat"), .combine = "rbind") %dopar% {
  seurat.obj = readRDS(paste0(dir.output, "/SeuratObject_", samples[j, "alias"], ".rds"))
  seurat.obj = Filter_ATAC(seurat = seurat.obj,
                           minPeakFrag = minPeakFrag,
                           maxPeakFrag = maxPeakFrag,
                           minPercPeak = minPercPeak,
                           maxNucleoSig = maxNucleoSig,
                           minTSS = minTSS)

  if(mito){
    seurat.obj = Add_MGATK(seurat = seurat.obj,
                           dir.data.sample = samples[j, "Directory"],
                           sample.ID = samples[j, "alias"],
                           MinCellVar = MinCellVar,
                           MinStrandCor = MinStrandCor,
                           MinVMR = MinVMR,
                           minMitoDepth = minMitoDepth)
  }

  saveRDS(seurat.obj, paste0(dir.output, "/SeuratObject_", samples[j, "alias"], ".rds"))
  seurat.obj@meta.data
}

if(cores > 1){
  parallel::stopCluster(cl)
}

if (mito){
  cols = setNames(c("#006600", "#006666", "#99CC99", "#CCCC33", "#CC9900"), c('pct_reads_in_peaks', 'peak_region_fragments', 'TSS.enrichment', 'nucleosome_signal', 'mtDNA_depth'))
} else{
  cols = setNames(c("#006600", "#006666", "#99CC99", "#CCCC33"), c('pct_reads_in_peaks', 'peak_region_fragments', 'TSS.enrichment', 'nucleosome_signal'))
  
}

pdf(paste0(dir.output, "/plots/vlnplot_qc_afterfiltering_ATAC.pdf"), width = (1+1.3*nrow(samples)), height = 4)
sapply(names(cols), function(x) {
  print(QCplot(df = meta, metric = x, color = cols[x], n.samples = nrow(samples), max.y = NULL, name.ID = name.ID))
})
dev.off()

if(mito){
  file.variants = list.files(dir.output, full.names = T)
  file.variants = grep("HighConf", file.variants, value = T)
  variants = unique(unlist(sapply(file.variants, function(x) read.table(x, header = T)$variant)))
  
  if(nb.cores > 1){
    cl <- parallel::makeCluster(nb.cores)
    doParallel::registerDoParallel(cl)
  } else {
    registerDoSEQ()
  }
  
  foreach(j = 1:nrow(samples), .packages = c("Signac", "Seurat"), .combine = "rbind") %dopar% {
    sample.ID = samples[j, "alias"]
    Calculate_Hetroplasmy(sample.ID = sample.ID, 
                          dir.output = dir.output,
                          variants = variants)
  }
  
  if(cores > 1){
    parallel::stopCluster(cl)
  }
}
