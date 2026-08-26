# =============================================================================
# 04_utility_tradeoff.R
# The core experiment: re-runs differential expression on each privatized
# (noised) count matrix, then measures how well the "significant genes"
# found in the private data overlap with the real ground-truth significant
# genes. This produces the privacy-utility tradeoff -- the central result
# of the project.
#
# Metrics used:
# - Jaccard index: |intersection| / |union| of the two significant gene sets
# - Precision: of the genes flagged significant under privacy, how many are
#   real (true positives)?
# - Recall: of the real significant genes, how many did privacy still let
#   us detect?
# =============================================================================

library(DESeq2)
library(SummarizedExperiment)

se               <- readRDS("data/processed/airway_se.rds")
sample_metadata  <- as.data.frame(colData(se))
ground_truth_sig <- readRDS("data/processed/ground_truth_sig_genes.rds")
epsilon_grid     <- readRDS("data/processed/epsilon_grid.rds")

run_de_on_matrix <- function(count_matrix, sample_metadata) {
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData   = sample_metadata,
    design    = ~ condition
  )
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds$condition <- relevel(dds$condition, ref = "Untreated")
  dds <- suppressMessages(DESeq(dds, quiet = TRUE))
  res <- results(dds, contrast = c("condition", "Treated", "Untreated"), alpha = 0.05)
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df$gene_id[
    !is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1
  ]
}

tradeoff_results <- data.frame(
  epsilon   = numeric(0),
  n_sig     = integer(0),
  jaccard   = numeric(0),
  precision = numeric(0),
  recall    = numeric(0),
  f1        = numeric(0)
)

for (eps in epsilon_grid) {
  fname <- paste0("data/processed/noised_matrices/counts_eps_", eps, ".rds")
  noised_counts <- readRDS(fname)

  cat("Running DE on epsilon =", eps, "...\n")
  sig_genes <- tryCatch(
    run_de_on_matrix(noised_counts, sample_metadata),
    error = function(e) character(0)
  )

  intersection <- length(intersect(sig_genes, ground_truth_sig))
  union_size   <- length(union(sig_genes, ground_truth_sig))
  jaccard      <- if (union_size == 0) NA else intersection / union_size
  precision    <- if (length(sig_genes) == 0) NA else intersection / length(sig_genes)
  recall       <- if (length(ground_truth_sig) == 0) NA else intersection / length(ground_truth_sig)
  f1           <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) {
    NA
  } else {
    2 * precision * recall / (precision + recall)
  }

  tradeoff_results <- rbind(tradeoff_results, data.frame(
    epsilon = eps, n_sig = length(sig_genes),
    jaccard = jaccard, precision = precision, recall = recall, f1 = f1
  ))
}

print(tradeoff_results)

dir.create("results", showWarnings = FALSE, recursive = TRUE)
write.csv(tradeoff_results, "results/privacy_utility_tradeoff.csv", row.names = FALSE)

cat("\nSaved to results/privacy_utility_tradeoff.csv. Now run 05_membership_inference.R.\n")
