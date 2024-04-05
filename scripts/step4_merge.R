#################################################################
#                                                               #
#                           STEP 4:                             #
#                       Merge samples                           #
#                                                               #
#################################################################
#################################################################
#                                                               #
#                  21.07.2023 - Ambre Giguelay                  #
#                                                               #
#################################################################

# -------- Load libraries -------- #
library(foreach)
library(ggplot2)
library(Signac)
library(Seurat)
library(harmony)

# -------- Load functions -------- #
message("Loading analysis functions")
source("scripts/signac_common.R")

# -------- Read parameters from config.yaml -------- #
dir.output = snakemake@params[['dir_sample']]
print(dir.output)
name.grp = "orig.ident"
print(name.grp)
GEX = snakemake@params[['GEX']]
print(GEX)
harmony = snakemake@params[['harmony']]
print(harmony)
CutOff_FTF = snakemake@params[['CutOff_FTF']]
print(CutOff_FTF)
NbDim = snakemake@params[['NbDim']]
print(NbDim)
RmCompo =snakemake@params[['RmCompo']]
print(RmCompo)
MinDistUMAP = snakemake@params[['MinDistUMAP']]
print(MinDistUMAP)
var.batch = snakemake@params[['var.batch']]
print(var.batch)
lambdaHarmony = snakemake@params[['lambdaHarmony']]
print(lambdaHarmony)
AlgoClustering = snakemake@params[['AlgoClustering']]
print(AlgoClustering)
ResolutionClustering = snakemake@params[['ResolutionClustering']]
print(ResolutionClustering)
batch_corr = snakemake@params[['batch_corr']]
print(batch_corr)

# -------- Read files, set variables -------- #
samples = read.table("config/samples.tsv", header = T)

col.palette = c("#333399", "#6666CC", "#9999FF", "#CCCCFF", "#e3b3e3", "#b37bb3", "#993366", "#82032d", "#bd2f00", "#e05b02", "#FF9900", "#e8ba02", "#b5b500", "#669900", "#52750c", "#009966", "cyan4", "turquoise3", "paleturquoise3", "powderblue", "lightskyblue3", "deepskyblue4", "lightsteelblue4", "#FFCC99", "#CC9966", "#996633", "#663300", "#CC9999", "#996666", "#FF9999", "#CC6666", "#e3acac", "#FFCCCC", "#FF99CC", "#FF6699", "#FF3399", "#CC3399", "black", "#0000CC", "#0066FF","#0099FF", "#65b5eb", "#95cbf0")

# -------- Run functions -------- #
files.seurat = list.files(dir.output)

if(!"SeuratObject_Merge.rds" %in% files.seurat){ #in case the merge has been done already
  seurat.list = setNames(sapply(unique(samples$condition), function(x) readRDS(paste0(dir.output, "/SeuratObjectBis_", x, ".rds"))), unique(samples$condition))
  
  print("Merging seurat objects")
  seurat = merge(x = seurat.list[[1]], y = seurat.list[-1])
  
  if(GEX){
    print("Computing GEX")
    DefaultAssay(seurat) = "peaks"
    # Compute gene accessibility
    gene.activities = GeneActivity(seurat)
    
    # Add to the Seurat object as a new assay and normalize
    seurat[['RNA']] = CreateAssayObject(counts = gene.activities)
    
    seurat = NormalizeData(
      object = seurat,
      assay = 'RNA',
      normalization.method = 'LogNormalize',
      scale.factor = median(seurat$nCount_RNA)
    )
  }
  message("Saving merged data")
  saveRDS(seurat, paste0(dir.output, "/SeuratObject_Merge.rds"))
  
} else {
  message("Reading merged object")
  seurat = readRDS(paste0(dir.output, "/SeuratObject_Merge.rds"))
}

if(harmony){
  message("Running harmony to correct batch effect")
  seurat@meta.data[[var.batch]] = as.factor(seurat@meta.data[[var.batch]])
  seurat = RunHarmony(object = seurat, group.by.vars = var.batch, reduction.use = 'lsi', assay.use = 'peaks', project.dim = FALSE, lambda = lambdaHarmony)
} else {
  message("Setting lsi reduction")
  reduction = "lsi"
}

message("Generating final seurat object")
seurat = SignacWorkflow(seurat = seurat,
                        dir.output = dir.output,
                        name.grp = name.grp, 
                        harmony = harmony,
                        CutOff_FTF = CutOff_FTF, 
                        NbDim = NbDim, 
                        RmCompo = RmCompo, 
                        Reduction = reduction,
                        MinDistUMAP = MinDistUMAP,
                        AlgoClustering = AlgoClustering,
                        ResolutionClustering = ResolutionClustering,
                        col.palette = col.palette)

saveRDS(seurat, paste0(dir.output, "/SeuratObject_Merge.rds"))
