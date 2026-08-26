# =============================================================================
# 05_membership_inference.R
# Simulates a simplified membership inference attack -- the real-world
# threat that differential privacy defends against. The question an
# attacker asks: "given a released statistic, was a specific person's data
# included in computing it?"
#
# THE ATTACK GAME (simplified, standard in DP teaching demos):
# For a given sample (the "target"), there are two possible worlds:
#   - "with":    the released statistic is the mean expression profile
#                INCLUDING the target
#   - "without": the released statistic is the mean expression profile
#                EXCLUDING the target
# The attacker is assumed to know both candidate true (noiseless) means
# (a standard, somewhat strong assumption in membership inference research,
# used here to demonstrate the mechanism clearly). They see ONE noised
# release and must guess which world produced it, by checking which
# candidate it's closer to.
#
# If privacy is weak (epsilon large / little noise), the released value is
# close to whichever true mean generated it, so the attacker guesses
# correctly almost every time. If privacy is strong (epsilon small / lots
# of noise), the two worlds become statistically indistinguishable, and
# attack accuracy drops toward 50% -- pure chance.
# =============================================================================

library(SummarizedExperiment)

se              <- readRDS("data/processed/airway_se.rds")
counts          <- assay(se)
clip_value      <- readRDS("data/processed/clip_value.rds")
epsilon_grid    <- readRDS("data/processed/epsilon_grid.rds")

# Restrict to expressed genes (same filter spirit as the DE analysis) so the
# attack is evaluated on genes that actually carry signal
keep <- rowSums(counts) >= 10
counts <- counts[keep, ]
cat("Running membership inference attack on", nrow(counts), "expressed genes.\n")

rlaplace <- function(n, scale) {
  u <- runif(n, -0.5, 0.5)
  -scale * sign(u) * log(1 - 2 * abs(u))
}

privatize_vector <- function(x, epsilon, clip_value) {
  if (is.infinite(epsilon)) return(x)
  clipped <- pmin(x, clip_value)
  clipped + rlaplace(length(clipped), scale = clip_value / epsilon)
}

n_trials_per_sample <- 100
attack_results <- data.frame(epsilon = numeric(0), accuracy = numeric(0))

set.seed(123)
for (eps in epsilon_grid) {
  correct <- 0
  total   <- 0

  for (i in seq_len(ncol(counts))) {
    true_mean_with    <- rowMeans(counts)
    true_mean_without <- rowMeans(counts[, -i, drop = FALSE])

    for (t in seq_len(n_trials_per_sample)) {
      true_label <- sample(c("with", "without"), 1)
      base <- if (true_label == "with") true_mean_with else true_mean_without
      released <- privatize_vector(base, eps, clip_value)

      d_with    <- sqrt(sum((released - true_mean_with)^2))
      d_without <- sqrt(sum((released - true_mean_without)^2))
      guess <- if (d_with < d_without) "with" else "without"

      correct <- correct + as.integer(guess == true_label)
      total   <- total + 1
    }
  }

  acc <- correct / total
  cat("epsilon =", eps, " -> attack accuracy:", round(acc, 3), "\n")
  attack_results <- rbind(attack_results, data.frame(epsilon = eps, accuracy = acc))
}

dir.create("results", showWarnings = FALSE, recursive = TRUE)
write.csv(attack_results, "results/membership_inference_results.csv", row.names = FALSE)

cat("\nSaved to results/membership_inference_results.csv.\n")
cat("Random-guess baseline is 0.5 -- accuracy near 0.5 means strong privacy protection.\n")
cat("Now run 06_visualization.R.\n")
