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
#                       Created by Ambre Giguelay               #
#    Adapted for snakemake by Coral Fustero-Torre               #
#                                                               #
#################################################################

# -------- Load libraries -------- #
library("foreach")
library("ggplot2")
library("Signac")
library("Seurat")

# -------- Load functions -------- #
message("Loading analysis functions")
source("scripts/signac_common.R")

# -------- Read parameters from config.yaml -------- #
dir.output = snakemake@params[['dir_sample']]
print(dir.output)
name.ID = "sample"
fc.resolution = snakemake@params[['fc_resolution']]
fc.k = snakemake@params[['fc_k']]
samples_tsv = snakemake@params[['samples_tsv']]
print(samples_tsv)
# -------- Read files, set variables -------- #
message("Reading metadata")
samples = read.table(samples_tsv, header = T)
message("Metadata loaded")
head(samples)
# -------- Run functions -------- #
cond = rev(unique(samples[, "condition"]))
print(cond)

foreach(grp = cond, .packages = c("Seurat", "Signac", "stringr")) %dopar%{
  samples.grp = samples[cond == grp,'sample']
  
  message(paste0("Calculating clonotypes for condition: ", cond)) 
          
  if(length(samples.grp) > 1){
    seurat.list = sapply(samples.grp, function(x) readRDS(paste0(dir.output, "/signac/SeuratObject_", x, ".s2.rds")))
    seurat = merge(
      x = seurat.list[[1]],
      y = seurat.list[names(seurat.list) != samples.grp[1]]
    )
    
  } else {
    seurat = readRDS(paste0(dir.output, "/signac/SeuratObject_", samples.grp, ".s2.rds"))
  }
  
  DefaultAssay(seurat) = "alleles"
  
  alleles = data.matrix(seurat[["alleles"]]@data)
  obj = getNN(t(sqrt(alleles)), k.param = fc.k)
  clusters = seuratSNN_cosineDistance(obj, resolution = fc.resolution) 
  cluster_name = str_pad(as.character(clusters), 3, pad = "0")
  
  seurat = AddMetaData(seurat, cluster_name, "Clonotype")
  
  saveRDS(seurat, paste0(dir.output, "/signac/SeuratObjectBis_", grp, ".rds"))
  
}

