# =============================================================================
# 01_load_data.R
# Loads the built-in "airway" RNA-seq dataset (airway smooth muscle cells,
# treated vs. untreated with dexamethasone, n=8). This dataset ships inside
# the airway Bioconductor package -- no download needed at all.
# =============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("airway", "DESeq2", "SummarizedExperiment"),
                      update = FALSE, ask = FALSE)

library(airway)
library(SummarizedExperiment)

data(airway)
se <- airway

# Clean up the condition column into something readable
se$condition <- factor(
  ifelse(se$dex == "trt", "Treated", "Untreated"),
  levels = c("Untreated", "Treated")
)

cat("Loaded airway dataset.\n")
cat("Samples:", ncol(se), " | Genes:", nrow(se), "\n")
print(table(se$condition))

# Save into this project's own data/processed/ folder
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
saveRDS(se, "data/processed/airway_se.rds")

cat("\nSaved to data/processed/airway_se.rds. Now run 02_ground_truth_de.R.\n")
