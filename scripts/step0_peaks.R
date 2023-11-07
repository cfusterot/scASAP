#################################################################
#                                                               #
#           Integrated quality control and Peak Calling         #
#                  Adapted from Ambre Giguelay's script         #
#                                by Coral Fustero-Torre         #
#                                                               #
#################################################################

# -------- Load libraries -------- #
suppressMessages(library("dplyr"))
suppressMessages(library("foreach"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("ggplot2"))
suppressMessages(library("Signac"))
suppressMessages(library("Seurat"))

# -------- Read parameters from config.yaml -------- #
message("Loading parameters")
min_peak_width = snakemake@params[["min_peak_width"]]
max_peak_width = snakemake@params[["max_peak_width"]]

# -------- Load functions -------- #
message("Loading analysis functions")
source("scripts/signac_common.R")

# -------- Set automatic parameters ------- #
message("Loading sample files")
outs = snakemake@params[["input_files"]]
dir_output = snakemake@params[["directory"]]
outdir = strsplit(dir_output, "/signac")[[1]]
samples_ID = unlist(sapply(strsplit(unlist(sapply(strsplit(outs, outdir), `[`, 2, simplify=FALSE)), "/cellranger_count/outs"), `[`, 1, simplify=FALSE))

# -------- Run functions -------- #
message("Generating common set of peaks for all samples:")
FindCommonPeaks_scATACseq(MaxPeakWidth = max_peak_width, MinPeakWidth = min_peak_width,
    dir.output = dir_output, dir.data.samples = outs)

