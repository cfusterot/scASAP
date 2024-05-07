library("GenomicAlignments")
library("SummarizedExperiment")
library("plyranges")
library("Rsamtools")
library("GenomeInfoDb")
library("GenomicRanges")
library("Biostrings")
library("BiocGenerics")
library("S4Vectors")
library("GenomicFeatures")
library("devtools")
if (!require("epiAneufinder")) {
  message("Installing epianeufinder software")
  install_github("colomemaria/epiAneufinder")
  library("epiAneufinder")
}

warnings()

message("Loading parameters")
fragments = snakemake@params["fragments"]
blacklist = snakemake@params["blacklist"]
windowSize = snakemake@params["windowSize"]
exclude = snakemake@params["exclude"]
output_dir = snakemake@params["output_dir"]

message("Loading genome")
genome = snakemake@params[["genome"]] 
if (genome == "hg38"){
  message("Loading BSgenome.Hsapiens.UCSC.hg38")
  suppressMessages(library("BSgenome.Hsapiens.UCSC.hg38"))
  genome = "BSgenome.Hsapiens.UCSC.hg38"
} else {
  message("Genome version not available.")
}

message("The following parameters will be used for running epianeufinder")
message(paste0("fragments: ", fragments))
message(paste0("blacklist: ", blacklist))
message(paste0("windowSize: ", windowSize))
message(paste0("exclude: ", exclude))
message(paste0("output_dir: ", output_dir))

epiAneufinder::epiAneufinder(input=fragments,
  outdir=output_dir,
  blacklist=blacklist,
  windowSize=1e5, 
  genome=genome,
  exclude=exclude, 
  reuse.existing=TRUE, 
  title_karyo="", 
  ncores=1, minFrags=20000)
