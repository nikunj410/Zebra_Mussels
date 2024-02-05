source("Bayesian_model/zebra/4_zebra_posterior_pred.R")
source("Bayesian_model/testing/4_synthetic_posterior_pred.R")

combined = plot_grid(First_detection_syn,First_detection_real, align = "v", ncol = 2)
plot(combined)
ggsave("Figures/Fig3_post_pred_check/poterior_predictive_checks.pdf",combined,  width = 6.5, height = 6, units = c("in"))
