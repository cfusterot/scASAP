library("Gotcha")

gotcha_JAK2_single = snakemake@params[["input_dir"]]
out_folder = snakemake@params[["output_dir"]]
# --- Analysis --- #
# 1. Merge mutation calls
print("Merging mutation calls")
MergeMutationCalling(out = out_folder)
