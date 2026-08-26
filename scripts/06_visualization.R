# =============================================================================
# 06_visualization.R
# Produces the two headline figures for this project:
# 1. The privacy-utility tradeoff curve (epsilon vs. F1/Jaccard)
# 2. The membership inference attack accuracy curve (epsilon vs. attack
#    success rate, with the 50% random-guess baseline marked)
# =============================================================================

library(ggplot2)

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# ---- Figure 1: Privacy-utility tradeoff ------------------------------------
tradeoff <- read.csv("results/privacy_utility_tradeoff.csv")

# Replace Inf with a large finite number just for plotting on a log axis
tradeoff$epsilon_plot <- ifelse(is.infinite(tradeoff$epsilon), 20, tradeoff$epsilon)

p_tradeoff <- ggplot(tradeoff, aes(x = epsilon_plot)) +
  geom_line(aes(y = f1, color = "F1 score"), linewidth = 1) +
  geom_point(aes(y = f1, color = "F1 score"), size = 2.5) +
  geom_line(aes(y = jaccard, color = "Jaccard index"), linewidth = 1, linetype = "dashed") +
  geom_point(aes(y = jaccard, color = "Jaccard index"), size = 2.5) +
  scale_x_continuous(trans = "log10") +
  scale_color_manual(values = c("F1 score" = "#2C6E9E", "Jaccard index" = "#B4483C")) +
  labs(
    title = "Privacy-Utility Tradeoff",
    subtitle = "How much differential expression accuracy survives at each privacy level",
    x = "Epsilon (privacy budget, log scale) \u2014 smaller = more private",
    y = "Score (vs. ground truth significant genes)",
    color = ""
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/privacy_utility_tradeoff.png", p_tradeoff, width = 7, height = 5, dpi = 150)

# ---- Figure 2: Membership inference attack accuracy ------------------------
attack <- read.csv("results/membership_inference_results.csv")
attack$epsilon_plot <- ifelse(is.infinite(attack$epsilon), 20, attack$epsilon)

p_attack <- ggplot(attack, aes(x = epsilon_plot, y = accuracy)) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "gray50") +
  annotate("text", x = min(attack$epsilon_plot), y = 0.52,
           label = "random-guess baseline", hjust = 0, size = 3, color = "gray40") +
  geom_line(color = "#8A4FA0", linewidth = 1) +
  geom_point(color = "#8A4FA0", size = 2.5) +
  scale_x_continuous(trans = "log10") +
  ylim(0.4, 1.0) +
  labs(
    title = "Membership Inference Attack Success Rate",
    subtitle = "Can an attacker tell if a sample was included in the released data?",
    x = "Epsilon (privacy budget, log scale) \u2014 smaller = more private",
    y = "Attacker accuracy"
  ) +
  theme_minimal()

ggsave("figures/membership_inference_attack.png", p_attack, width = 7, height = 5, dpi = 150)

cat("Figures saved: figures/privacy_utility_tradeoff.png, figures/membership_inference_attack.png\n")
cat("\n=== Project complete ===\n")
