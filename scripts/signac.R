#################################################################
#                                                               #
#           Integrated quality control and Peak Calling         #
#                  Adapted from Ambre Giguelay's script         #
#                                by Coral Fustero-Torre         #
#                                                               #
#################################################################

# -------- Load libraries -------- #
suppressMessages(library("Biostrings"))
suppressMessages(library("BSgenome.Hsapiens.UCSC.hg38")) #BSgenome.Mmusculus.UCSC.mm10
suppressMessages(library("dplyr"))
suppressMessages(library("EnsDb.Hsapiens.v86")) #for mm10
suppressMessages(library("foreach"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("ggplot2"))
suppressMessages(library("harmony"))
suppressMessages(library("Signac"))
suppressMessages(library("Seurat"))
suppressMessages(library("data.table"))

# -------- Load functions -------- #
message("Loading analysis functions")
source("scripts/common_signac.R")

# -------- Folder configuration  -------- #
#folders = c("qc", "mito", "visualisation")

# -------- Read parameters from config.yaml -------- #
message("Loading parameters")
min_peak_width = snakemake@params[["min_peak_width"]]
max_peak_width = snakemake@params[["max_peak_width"]]
min_counts = snakemake@params[["min_counts"]]
min_peak_fragment = snakemake@params[["min_peak_fragment"]]
max_peak_fragment = snakemake@params[["max_peak_fragment"]]
max_nucleosome_signal = snakemake@params[["max_nucleosome_signal"]]
min_TSS = snakemake@params[["min_TSS"]]
min_percentage_peaks = snakemake@params[["min_percentage_peaks"]]
# Mito data
min_depth = snakemake@params[["min_depth"]]
min_cell_var = snakemake@params[["min_cell_var"]]
min_VMR = snakemake@params[["min_VMR"]]
min_strand_cor = snakemake@params[["min_strand_cor"]]

# -------- Set automatic parameters ------- #
message("Loading sample files")

outs = snakemake@input[["outs"]]
mgatk = snakemake@input[["mgatk"]]
amulet = snakemake@input[["amulet"]]
dir_output = snakemake@output[["directory"]]
outdir = strsplit(dir_output, "/signac")[[1]]
samples.ID = unlist(sapply(strsplit(unlist(sapply(strsplit(outs, outdir), `[`, 2, simplify=FALSE)), "/cellranger_count"), `[`, 1, simplify=FALSE))

# -------- Run functions -------- #
message("START:")
message("Step1.- Computing set of common peaks")
mylists = FindCommonPeaks_scATACseq(samples.ID = samples.ID,
                                    granges.ens = EnsDb.Hsapiens.v86,
                                    genome = "hg38",
                                    dir.data.samples = setNames(outs, samples.ID),
                                    dir.output =dir_output,
                                    dir.mito = setNames(mgatk, samples.ID), MaxPeakWidth = max_peak_width, # default: 10000
                                    MinPeakWidth = min_peak_width, # default: 20
                                    MinCounts = min_counts, #default: 500
                                    analysis.ID = "Batch1",
                                    name.ID = "Dataset",
                                    saveRDS = T, 
                                    name.rds = "SeuratObjectRaw", 
                                    mito = T)
message("Step2.- Adding AMULET output to Signac object")
mylists[["Doublets"]] = MakeDoubletList_scATACseq(samples.ID = samples.ID ,dir.doublets = amulet)

message("Step2.- Saving filtered results")
mylists[["Seurat"]] = lapply(names(mylists[["Seurat"]]), function(name) {FilterCells_scATACseq(seurat.obj = mylists[["Seurat"]][[name]],
                                                                                 samples.ID = name, 
                                                                                 analysis.ID = "Batch1",
                                                                                 dir.output = dir_output,
                                                                                 saveRDS = T, 
                                                                                 name.rds = "SeuratObjectFilter", 
                                                                                 doublets = mylists[["Doublets"]][[name]], 
                                                                                 minPeakFrag = min_peak_fragment, #default: 1000
                                                                                 maxPeakFrag = max_peak_fragment, #default: 30000
                                                                                 minPercPeak = min_percentage_peaks, #default: 30
                                                                                 maxNucleoSig = max_nucleosome_signal, #default: 2
                                                                                 minTSS = min_TSS, #default: 2
                                                                                 minMitoDepth = min_depth, #default: 10
                                                                                 mito = T)})


message("Step3.- Identifying mitochondrial variants")
names(mylists[["Seurat"]]) = names(mylists[["MGATK"]])
mylists[["Variant"]] = lapply(names(mylists[["Seurat"]]), function(name) {IdentVariant_scATACseq(seurat.obj = mylists[["Seurat"]][[name]],
                                                                                  mgatk = mylists[["MGATK"]][[name]],
                                                                                  samples.ID = name,
                                                                                  analysis.ID = "Batch1",
                                                                                  dir.output = dir_output,
                                                                                  saveVariant = T,
                                                                                  MinCellVar = min_cell_var #default: 5
                                                                                  )})
mylists[["MGATK"]] = NULL
names(mylists[["Variant"]]) = names(mylists[["Seurat"]])

message("Step4.- Computing Clonotypes")
mylists[["Seurat"]] = lapply(names(mylists[["Seurat"]]), function(name) {FindClonotypes_scATACseq(seurat.obj = mylists[["Seurat"]][[name]], 
                                                                                      variable.sites = mylists[["Variant"]][[name]], 
                                                                                      samples.ID = name, 
                                                                                      dir.output = dir_output, 
                                                                                      saveRDS = T, 
                                                                                      name.rds = "SeuratObjectClonotypes",
                                                                                      analysis.ID = "Batch1",
                                                                                      MinStrandCor = min_strand_cor, #default: 0.65, 
                                                                                      MinCellVar = min_cell_var, #default: 5
                                                                                      MinVMR = min_VMR #default: 0.01
                                                                                      )})
names(mylists[["Seurat"]]) = names(mylists[["Variant"]])
mylists[["Variant"]] = NULL

message("Step5.- Merging final result")
seurat.all = Merge_scATACseq(mylists[["Seurat"]], 
                             dir.output = dir_output, 
                             saveRDS = T, 
                             name.rds = "SeuratCombined",
                             analysis.ID = "Batch1",
                             name.ID = "Dataset", 
                             mito = T)

message("Step5.- Calculating Gene activity")
seurat.all = GeneActivity_scATACseq(seurat.obj = seurat.all,
                                    dir.output = dir_output,
                                    saveRDS = T, 
                                    name.rds = "SeuratCombined_GeneScore",
                                    analysis.ID = "Batch1")
