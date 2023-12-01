#################################################################
#                                                               #
#           Integrated quality control and Peak Calling         #
#                  Adapted from Ambre Giguelay's script         #
#                                by Coral Fustero-Torre         #
#                                                               #
#################################################################

# -------- Load libraries -------- #
message("Loading libraries")
suppressMessages(library("dplyr"))
suppressMessages(library("foreach"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("ggplot2"))
suppressMessages(library("Signac"))
suppressMessages(library("Seurat"))

message("Loading genome")
granges_ens = snakemake@params[["granges_ens"]]
genome = snakemake@params[["genome"]] 
if (grepl("Hsapiens", granges_ens, fixed = TRUE)){
  message("Loading EnsDb.Hsapiens.v86")
  suppressMessages(library("EnsDb.Hsapiens.v86"))
  gene_annotation = GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  gene_annotation = renameSeqlevels(gene_annotation, mapSeqlevels(seqlevels(gene_annotation), "UCSC"))
  genome(gene_annotation) = genome
} else if (grepl("Mmusculus", granges_ens, fixed = TRUE)) {
  message("Loading EnsDB.Mmusculus.v79")
  suppressMessages(library("EnsDb.Mmusculus.v79"))
  gene_annotation = GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  gene_annotation = renameSeqlevels(gene_annotation, mapSeqlevels(seqlevels(gene_annotation), "UCSC"))
} else {
    message("Species not found.")
}

# -------- Read parameters from config.yaml -------- #
message("Loading parameters")
min_peak_width = snakemake@params[["min_peak_width"]]
max_peak_width = snakemake@params[["max_peak_width"]]
min_counts = snakemake@params[["min_counts"]]
integration = snakemake@params[["integration"]] 
nb_cores = snakemake@params[["nb_cores"]]
genome = snakemake@params[["genome"]]

# -------- Load functions -------- #
message("Loading analysis functions")
source("scripts/signac_common.R")

# -------- Set automatic parameters ------- #
message("Loading sample files")
mgatk = snakemake@input[["mgatk"]]
amulet = snakemake@input[["amulet"]]
dir_integration = snakemake@params[["dir_integration"]]
dir_sample = snakemake@params[["dir_sample"]]
sample_ID = snakemake@params[["sample_ID"]]
sample_ID = strsplit(sample_ID, '[{}]')[[1]][2]

# -------- Run functions -------- #
message("START:")
message(paste0("Creating ATAC object for sample: ", sample_ID))
meta <- CreateATACobject(sample.ID = sample_ID,
                   dir.data.sample = dir_sample,
                   MinCounts = min_counts,
                   dir.output = dir_integration,
                   gannotation = gene_annotation,
                   name.ID = sample_ID)

cols = setNames(c("#006600", "#006666", "#99CC99", "#CCCC33"), c('pct_reads_in_peaks', 'peak_region_fragments', 'TSS.enrichment', 'nucleosome_signal'))

# -- Create plots directory -- #
dir.create(paste0(dir_integration, "plots/"))
message("Generating basic QC plots")
pdf(paste0(dir_integration, "/plots/",sample_ID,"_vlnplot_qc_beforeFiltering_ATAC.pdf"), width = 8, height = 4)
sapply(names(cols), function(x) {
  print(x)
  QCplot(df = meta, metric = x, color = cols[x], n.samples = 1, max.y = NULL, name.ID = sample_ID)
})
dev.off()

