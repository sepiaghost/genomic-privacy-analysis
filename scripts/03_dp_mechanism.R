# =============================================================================
# 03_dp_mechanism.R
# Implements the Laplace mechanism (the standard differential privacy
# technique) and applies it to the gene expression count matrix at several
# privacy levels (epsilon). Produces one "privatized" count matrix per
# epsilon value.
#
# HOW IT WORKS (the short version):
# - Differential privacy adds carefully calibrated random noise to data so
#   that no single individual's presence or absence can be reliably
#   detected from the output, while still preserving overall statistical
#   patterns.
# - "Epsilon" controls the privacy budget: SMALLER epsilon = MORE noise =
#   MORE privacy = LESS accuracy. Epsilon = Inf means no privacy protection
#   at all (the original, unmodified data).
# - The Laplace mechanism adds noise drawn from a Laplace distribution,
#   scaled by (sensitivity / epsilon). "Sensitivity" is the maximum amount
#   one individual's data could possibly change the output. We bound this
#   by CLIPPING every count to a maximum value first -- otherwise a single
#   extremely highly-expressed gene in one sample could have unbounded
#   influence, which would break the privacy guarantee.
#
# HONEST SCOPING NOTE: applying independent per-gene, per-sample noise like
# this is illustrative of the core mechanism, not a formally composed
# privacy guarantee for releasing an entire count matrix (that requires
# careful privacy budget accounting across every gene and sample queried --
# see "composition" in the differential privacy literature). This project
# demonstrates the mechanism and its utility tradeoff, not a
# production-grade private release system.
# =============================================================================

library(SummarizedExperiment)

se <- readRDS("data/processed/airway_se.rds")
counts <- assay(se)

# ---- Laplace distribution sampler (base R has no built-in rlaplace) -------
rlaplace <- function(n, scale) {
  u <- runif(n, -0.5, 0.5)
  -scale * sign(u) * log(1 - 2 * abs(u))
}

# ---- The Laplace mechanism --------------------------------------------------
laplace_mechanism <- function(count_matrix, epsilon, clip_value) {
  if (is.infinite(epsilon)) {
    return(count_matrix)  # epsilon = Inf means "no privacy" -- return as-is
  }
  clipped <- pmin(count_matrix, clip_value)
  noise_scale <- clip_value / epsilon
  noise <- matrix(
    rlaplace(length(clipped), scale = noise_scale),
    nrow = nrow(clipped), ncol = ncol(clipped)
  )
  noised <- clipped + noise
  # Post-processing: counts can't be negative or fractional. This is safe
  # under DP's post-processing property -- it doesn't weaken the guarantee.
  noised <- pmax(round(noised), 0)
  dimnames(noised) <- dimnames(count_matrix)
  noised
}

# ---- Choose a clipping bound -----------------------------------------------
# Use the 95th percentile of all counts as a reasonable clip value -- bounds
# sensitivity without destroying the signal from most genes.
clip_value <- quantile(counts, 0.95)
cat("Clipping bound (95th percentile of counts):", clip_value, "\n")

# ---- Apply the mechanism across a range of privacy levels ------------------
epsilon_grid <- c(0.1, 0.5, 1, 2, 5, 10, Inf)

dir.create("data/processed/noised_matrices", showWarnings = FALSE, recursive = TRUE)

set.seed(42)  # reproducible noise draws
for (eps in epsilon_grid) {
  noised <- laplace_mechanism(counts, epsilon = eps, clip_value = clip_value)
  fname <- paste0("data/processed/noised_matrices/counts_eps_", eps, ".rds")
  saveRDS(noised, fname)
  cat("epsilon =", eps, " -> saved to", fname, "\n")
}

saveRDS(clip_value, "data/processed/clip_value.rds")
saveRDS(epsilon_grid, "data/processed/epsilon_grid.rds")

cat("\nAll privatized matrices generated. Now run 04_utility_tradeoff.R.\n")
