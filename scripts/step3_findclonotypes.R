#################################################################
#                                                               #
#                Profiling of PolG mice tissues                 #          
#                                                               #
#                           STEP 3:                             #
#                     Define clonotypes                         #
#                                                               #
#################################################################
#################################################################
#                                                               #
#                  21.07.2023 - Ambre Giguelay                  #
#                                                               #
#################################################################
rm(list = ls())

# -------- Load libraries -------- #
library("foreach")
library("ggplot2")
library("Signac")
library("Seurat")

# -------- Load functions -------- #
message("Loading analysis functions")
source("scripts/signac_common.R")

# -------- Read parameters from config.yaml -------- #
dir.output = snakemake@output[["directory"]]
nb.cores = snakemake@params['cores']
name.ID = "sample"
fc.resolution = snakemake@params['fc_resolution']
fc.k = snakemake@params['fc_k']

# -------- Read files, set variables -------- #
samples = read.table("config/samples.tsv", header = T)

# -------- Run functions -------- #
if(nb.cores > 1){
  cl <- parallel::makeCluster(nb.cores)
  doParallel::registerDoParallel(cl)
} else {
  registerDoSEQ()
}

options(future.globals.maxSize = 8000 * 1024^2) #To avoid error message with large datasets

foreach(grp = rev(unique(samples$condition)), .packages = c("Seurat", "Signac", "stringr")) %dopar%{
  samples.grp = samples[samples$condition == grp,]$alias
  
  if(length(samples.grp) > 1){
    seurat.list = sapply(samples.grp, function(x) readRDS(paste0(dir.output, "/SeuratObject_", x, ".rds")))
    seurat = merge(
      x = seurat.list[[1]],
      y = seurat.list[names(seurat.list) != samples.grp[1]]
    )
    
  } else {
    seurat = readRDS(paste0(dir.output, "/SeuratObject_", samples.grp, ".rds"))
  }
  
  DefaultAssay(seurat) = "alleles"
  
  alleles = data.matrix(seurat[["alleles"]]@data)
  obj = getNN(t(sqrt(alleles)), k.param = fc.k)
  clusters = seuratSNN_cosineDistance(obj, resolution = fc.resolution) 
  cluster_name = str_pad(as.character(clusters), 3, pad = "0")
  
  seurat = AddMetaData(seurat, cluster_name, "Clonotype")
  
  saveRDS(seurat, paste0(dir.output, "/SeuratObjectBis_", grp, ".rds"))
  
}

if(nb.cores > 1){
  parallel::stopCluster(cl)
}