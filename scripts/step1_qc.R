#################################################################
#                                                               #
#           Integrated quality control and Peak Calling         #
#                  Adapted from Ambre Giguelay's script         #
#                                by Coral Fustero-Torre         #
#                                                               #
#################################################################

# -------- Load libraries -------- #
granges_ens = snakemake@params[["granges_ens"]]
if (grepl("Hsapiens", ranges, fixed = TRUE){
  message("Loading EnsDb.Hsapiens.v86")
  suppressMessages(library("EnsDb.Hsapiens.v86"))
} else if (grepl("Mmusculus", ranges, fixed = TRUE)) {
  message("Loading EnsDB.Mmusculus.v79")
  suppressMessages(library("EnsDb.Mmusculus.v79"))
} else {
  message("Species not found.")
}
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
min_counts = snakemake@params[["min_counts"]]
integration = snakemake@params[["integration"]] 
nb_cores = snakemake@params[["nb_cores"]]
granges_ens = snakemake@params[["granges_ens"]] 
genome = snakemake@params[["genome"]]

# -------- Load functions -------- #
message("Loading analysis functions")
source("scripts/signac_common.R")

# -------- Set automatic parameters ------- #
message("Loading sample files")
outs = snakemake@input[["outs"]]
mgatk = snakemake@input[["mgatk"]]
amulet = snakemake@input[["amulet"]]
dir_output = snakemake@output[["directory"]]
outdir = strsplit(dir_output, "/signac")[[1]]
samples_ID = unlist(sapply(strsplit(unlist(sapply(strsplit(outs, outdir), `[`, 2, simplify=FALSE)), "/cellranger_count/outs"), `[`, 1, simplify=FALSE))

# -------- Run functions -------- #
message("START:")
gene_annotation = GetGRangesFromEnsDb(ensdb = granges_ens)
gene_annotation = renameSeqlevels(gene_annotation, mapSeqlevels(seqlevels(gene_annotation), "UCSC"))
genome(gene_annotation) = genome

message(paste0("Creating ATAC object for sample: ", samples_ID))
CreateATACobject(sample.ID = samples_ID,
                   dir.data.sample = outs,
                   MinCounts = min_counts,
                   dir.output = dir_output,
                   gannotation = gene_annotation,
                   name.ID = "sample")

cols = setNames(c("#006600", "#006666", "#99CC99", "#CCCC33"), c('pct_reads_in_peaks', 'peak_region_fragments', 'TSS.enrichment', 'nucleosome_signal'))
dir.create(paste0(outdir, "plots/"))

message("Generating basic QC plots")
pdf(paste0(dir.output, "/plots/vlnplot_qc_beforeFiltering_ATAC.pdf"), width = (1+1.3*nrow(samples)), height = 4)
QCplot(df = meta, metric = x, color = cols[x], n.samples = nrow(samples), max.y = NULL, name.ID = name.ID)
dev.off()

