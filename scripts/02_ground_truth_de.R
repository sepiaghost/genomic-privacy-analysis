# =============================================================================
# 02_ground_truth_de.R
# Establishes the "ground truth" -- the real, non-private differential
# expression result on the original data. Every privacy-noised version
# later gets compared back against this to measure utility loss.
# =============================================================================

library(DESeq2)
library(SummarizedExperiment)

se <- readRDS("data/processed/airway_se.rds")

counts <- assay(se)
sample_metadata <- as.data.frame(colData(se))

# ---- Standard DESeq2 differential expression on the real data -------------
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = sample_metadata,
  design    = ~ condition
)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds$condition <- relevel(dds$condition, ref = "Untreated")
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "Treated", "Untreated"), alpha = 0.05)
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# ---- Define the "ground truth" significant gene set ------------------------
ground_truth_sig <- res_df$gene_id[
  !is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1
]

cat("Ground truth: ", length(ground_truth_sig),
    " significant genes (padj<0.05, |log2FC|>1) out of ", nrow(res_df),
    " tested.\n", sep = "")

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

write.csv(res_df, "results/ground_truth_de_results.csv", row.names = FALSE)
saveRDS(ground_truth_sig, "data/processed/ground_truth_sig_genes.rds")
saveRDS(dds, "data/processed/ground_truth_dds.rds")

cat("\nSaved. Now run 03_dp_mechanism.R.\n")
