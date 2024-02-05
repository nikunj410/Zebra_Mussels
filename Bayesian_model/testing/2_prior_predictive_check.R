# Running 4 chains for prior sampling takes 20 secs to run on Macbook pro 2022 (Apple M1 Max) 

library(cmdstanr)
library(readr)
library(stringr)
library(rethinking)

histospark <- function(x, width = 11) {
  sparks <- c("\u2581", "\u2582", "\u2583","\u2584", "\u2585","\u2586", "\u2587")
  bins <- graphics::hist(x, breaks = seq(min(x),max(x),length.out=width), plot = FALSE)
  factor <- cut(
    bins$counts / max(bins$counts),
    breaks = seq(0, 1, length = length(sparks) + 1),
    labels = sparks,
    include.lowest = TRUE
  )
  paste0(factor, collapse = "")
}
true_index <- function(x, true, width = 11) {
  sparks <- c("\u2581", "\u2582", "\u2583", "\u2585","\u2586", "\u2587")
  bins <- graphics::hist(x, breaks = seq(min(x),max(x),length.out=width), plot = FALSE)
  return(tail(which(bins$breaks<true),1))
}

options(mc.cores = parallel::detectCores())

Adj = read_csv("Bayesian_model/processing_raw_data/data_files/adj/Adj_hex.csv",show_col_types = FALSE)
survey = as.matrix(read_csv("Bayesian_model/testing/data_files/synthetic_survey.csv", show_col_types = FALSE))
sigma =  read_csv("Bayesian_model/processing_raw_data/data_files/life_history_prior/sigma.csv", show_col_types = FALSE)
f =  read_csv("Bayesian_model/processing_raw_data/data_files/life_history_prior/f.csv", show_col_types = FALSE)
DistM = read_csv("Bayesian_model/processing_raw_data/data_files/DistM/DistM.csv",show_col_types = FALSE)

A <- as.matrix(Adj/sum(Adj))
rownames(A)=colnames(A)

dat <- list(n = ncol(survey), Time = nrow(survey), A = A, survey = survey, Dist = DistM, 
            n_port = ncol(survey), index_port = 1:ncol(survey),
            nonzero = (length(which(A!=0))+ncol(survey)),n_ages = 4,
            sigma_mu = (sigma$Max+sigma$Min)/2, sigma_sd = (sigma$Max-sigma$Min)/5, N_star = 1e-9,
            f_mu = (f$Max+f$Min)/4, f_sd = (f$Max-f$Min)/10,
            prior_sampling = as.integer(TRUE), posterior_checks = as.integer(FALSE))

file <- file.path("Bayesian_model/Stan_model/", "zebra_fit.stan") # with beta distribution
mod <- cmdstan_model(file, stanc_options = list("O1"))
path = c("Bayesian_model/testing/data_files/synthetic_prior_sampling_stan.rds")

fit <- mod$sample(
  data = dat,
  chains = 4,
  seed = 122,
  parallel_chains = 8,
  show_messages = T,
  show_exceptions = F,
  refresh = 100,
  iter_warmup = 300,
  iter_sampling = 700
)

parameter_draw = as.data.frame(fit$draws( format="df", variables = c("r", "Deff", "intro_year", "N_knot","p_Kd","p_mus", "p_sampling")))
parameter_draw = parameter_draw[,-((ncol(parameter_draw)-2):ncol(parameter_draw))]
parameter_summary = precis( parameter_draw, prob=0.95, depth=2)

parameter_summary$variable = rownames(parameter_summary)
rownames(parameter_summary) <- NULL
parameter_summary = parameter_summary[, c(ncol(parameter_summary),1,3:(ncol(parameter_summary)-1))]
parameter_summary[,2:(ncol(parameter_summary)-1)] = format(round(parameter_summary[,2:(ncol(parameter_summary)-1)], 3), nsmall = 3)
parameter_summary$histogram = sapply(1:ncol(parameter_draw),function(x) histospark(parameter_draw[,x]))
print(parameter_summary)
write_csv(parameter_summary, "Bayesian_model/testing/data_files/synthetic_prior_pred_checks.csv")




