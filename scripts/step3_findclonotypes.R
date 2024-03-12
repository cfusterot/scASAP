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
message(paste0("Setting output directory to: ", dir.output))


name.ID = snakemake@params[['sample_ID']]
fc.resolution = snakemake@params[['fc_resolution']]
fc.k = snakemake@params[['fc_k']]
message("Setting analysis parameters")

samples_tsv = snakemake@params[['samples_tsv']]
message(paste0("Setting samples.tsv directory to: ", samples_tsv))

# -------- Read files, set variables -------- #
message("Loading metadata table")
samples = read.table(samples_tsv, header = T, row.names = 1)
message(samples)

# -------- Run functions -------- #
cond = unique(samples[name.ID, "condition"])
message(paste0("Calculating clonotypes for condition: ", cond)) 
          
if(length(name.ID) > 1){
  seurat.list = sapply(seq_along(name.ID), function(x) readRDS(paste0(dir.output[x], "/signac/SeuratObject_", name.ID[x], ".s2.rds")))
  message("Merging samples from the same condition")
  seurat = merge(x = seurat.list[[1]], y = seurat.list[names(seurat.list) != name.ID[1]])
} else {
  seurat = readRDS(paste0(dir.output, "/signac/SeuratObject_", name.ID, ".s2.rds"))
}
  
DefaultAssay(seurat) = "alleles"
  
alleles = data.matrix(seurat[["alleles"]]@data)
obj = getNN(t(sqrt(alleles)), k.param = fc.k)
clusters = seuratSNN_cosineDistance(obj, resolution = fc.resolution) 
cluster_name = str_pad(as.character(clusters), 3, pad = "0")
  
seurat = AddMetaData(seurat, cluster_name, "Clonotype")

saveRDS(seurat, paste0(dir.output, "/signac/SeuratObjectBis_", cond, ".rds"))

