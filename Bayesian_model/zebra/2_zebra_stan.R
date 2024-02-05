# Running 4 MCMC chains take ~23 hours to run on Macbook pro 2022 (Apple M1 Max) if ones uses all the ports as prior for introduction (# All 160 ports)
# Running 4 MCMC chains take ~2.25 hours to run on Macbook pro 2022 (Apple M1 Max) if one uses ports around lake Erie as prior for introduction (# 14 ports)

# Running 4 multi-pathfinder chains take ~400 secs to run on Macbook pro 2022 (Apple M1 Max) if ones uses all the ports as prior for introduction (# All 160 ports)
# Running 4 multi-pathfinder chains take ~40 secs to run on Macbook pro 2022 (Apple M1 Max) if one uses ports around lake Erie as prior for introduction (# 14 ports)

# The model takes ~19 secs to compile on Macbook pro 2022 (Apple M1 Max)

library(cmdstanr)
library(readr)
library(rtern)

options(mc.cores = parallel::detectCores())

Adj = read_csv("Bayesian_model/processing_raw_data/data_files/adj/Adj_hex.csv",show_col_types = FALSE)
survey = as.matrix(read_csv("Bayesian_model/zebra/data_files/zebra_survey_matrix.csv", show_col_types = FALSE))
sigma =  read_csv("Bayesian_model/processing_raw_data/data_files/life_history_prior/sigma.csv", show_col_types = FALSE)
f =  read_csv("Bayesian_model/processing_raw_data/data_files/life_history_prior/f.csv", show_col_types = FALSE)
intro_port_prior = read_csv("Bayesian_model/processing_raw_data/data_files/intro_port_prior/intro_port_prior_hexID.csv", show_col_types = FALSE)
DistM = read_csv("Bayesian_model/processing_raw_data/data_files/DistM/DistM.csv", show_col_types = FALSE)

A <- as.matrix(Adj/sum(Adj))
rownames(A)=colnames(A)
intro_port_index_prior = match(intro_port_prior$prior_intro_port,as.integer(colnames(A)))
Erie_as_prior = 1

dat <- list(n = ncol(survey), Time = nrow(survey), A = A, survey = survey, Dist = DistM,
            n_port = Erie_as_prior == 1 ? length(intro_port_index_prior) : ncol(survey),
            index_port = Erie_as_prior == 1 ? intro_port_index_prior : (1:ncol(survey)),
            nonzero = (length(which(A!=0))+ncol(survey)),n_ages = 4,
            sigma_mu = (sigma$Max+sigma$Min)/2, sigma_sd = (sigma$Max-sigma$Min)/5,
            f_mu = (f$Max+f$Min)/4, f_sd = (f$Max-f$Min)/10, posterior_checks = as.integer(TRUE),
            prior_sampling = as.integer(FALSE))

file <- file.path("Bayesian_model/Stan_model/", "zebra_fit.stan")
mod <- cmdstan_model(file, stanc_options = list("O1"), cpp_options = list(stan_threads = TRUE))

fit_pf_multi <- mod$pathfinder(data = dat, num_paths = 4, single_path_draws = 40, history_size=50, max_lbfgs_iters=100, num_threads = 10)
print(fit_pf_multi$summary(c("r", "Deff", "intro_year", "p_Kd","p_mus", "log_No", "h", "sigma", "f", "K_hyper")))
mcmc_hist(fit_pf_multi$draws("intro_port"))
initial_opt = as.vector(colMeans(as.data.frame(fit_pf_multi$draws( format="df", variables = c("lambda","log_No","h", "deff", "beta_mu_obs", "K_hyper", "sigma_std", "f_std", "beta_obs")))))
init_fun <- function() list(lambda = initial_opt[1], log_No = initial_opt[2], h=initial_opt[3], deff = initial_opt[4], beta_mu_obs = initial_opt[5],
                            K_hyper = initial_opt[6:8], sigma_std = initial_opt[9:13], f_std = initial_opt[14:16], beta_obs = initial_opt[17:176])


path = path = c("Bayesian_model/zebra/data_files/zebra_stan.rds")

fit <- mod$sample(
  data = dat,
  chains = 4,
  seed = 123,
  parallel_chains = 8,
  show_messages = T,
  show_exceptions = F,
  init = init_fun,
  refresh = 100,
  iter_warmup = 300,
  threads_per_chain = 1,
  iter_sampling = 700
)
print(fit$summary(c("r", "Deff", "intro_year", "p_Kd","p_mus", "log_No", "h", "sigma", "f", "K_hyper")))
mcmc_hist(fit$draws("intro_port"))
fit$save_object(file = path)


