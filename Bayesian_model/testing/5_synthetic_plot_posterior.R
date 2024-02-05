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

mc.cores = parallel::detectCores()
source("Bayesian_model/testing/0_synthetic_parameters.R")
fit = readRDS("Bayesian_model/testing/data_files/synthetic_stan.rds")
HexID = read_csv("Bayesian_model/processing_raw_data/data_files/port_list/portsHexID.csv", show_col_types = FALSE)
Adj = read_csv("Bayesian_model/processing_raw_data/data_files/adj/Adj_hex.csv",show_col_types = FALSE)
sigma =  read_csv("Bayesian_model/processing_raw_data/data_files/life_history_prior/sigma.csv", show_col_types = FALSE,col_names = TRUE)
f =  read_csv("Bayesian_model/processing_raw_data/data_files/life_history_prior/f.csv", show_col_types = FALSE,col_names = TRUE)
summary = fit$summary(c("r", "Deff", "intro_year", "p_Kd","p_mus", "log_No", "h", "sigma", "f","K_hyper","p_sampling"))

intro_draw = fit$draws(variables = "intro_port", format="df")
intro_table = as.data.frame(table(as.integer(colnames(Adj)[intro_draw$intro_port])))
intro_table$Freq = as.numeric(format(round(intro_table$Freq/sum(intro_table$Freq), 3), nsmall = 3))*100
intro_table = intro_table[rev(order(intro_table$Freq)),]
colnames(intro_table)  = c("hex_id","percent")
intro_table$lon = HexID$hex_lon[match(intro_table$hex_id,HexID$hex_id)]
intro_table$lat = HexID$hex_lat[match(intro_table$hex_id,HexID$hex_id)]
print(intro_table)
write_csv(as.data.frame(intro_table),"Bayesian_model/testing/data_files/synthetic_intro_port.csv")

true_values = c(r, Deff, IntroYear, P_Kd, P_mus, log10(N_knot), H, (sigma$Max+sigma$Min)/2, (f$Max+f$Min)/4, Eta, Rho, Sigma, p_sampling)
parameter_draw = as.data.frame(fit$draws( format="df", variables = c("r", "Deff", "intro_year", "p_Kd","p_mus", "log_No", "h", "sigma", "f","K_hyper","p_sampling")))
parameter_draw = parameter_draw[,-((ncol(parameter_draw)-2):ncol(parameter_draw))]

parameter_draw = rbind(parameter_draw,true_values)
parameter_summary = cbind(precis( parameter_draw, prob=0.95, depth=2), summary[,c(8,9)])

parameter_summary$true_value = true_values
pass_fit_test = true_values>parameter_summary$`2.5%` & true_values<parameter_summary$`97.5%`
parameter_summary$variable = rownames(parameter_summary)
rownames(parameter_summary) <- NULL
parameter_summary = parameter_summary[, c(9,8,1,3,4,6,7,5)]
parameter_summary[,2:(ncol(parameter_summary)-1)] = format(round(parameter_summary[,2:(ncol(parameter_summary)-1)], 2), nsmall = 2)
parameter_summary$histogram = sapply(1:ncol(parameter_draw),function(x) histospark(parameter_draw[,x]))
parameter_summary$true_index = sapply(1:ncol(parameter_draw),function(x) true_index(parameter_draw[,x],true_values[x]))
parameter_summary$pass_fit_test = sapply(pass_fit_test, function(x){if (x){c("Y")} else {c("N")}})
parameter_summary$chains_convergence = sapply(as.numeric(parameter_summary$rhat)>1.05, function(x){if (x){c("N")} else {c("Y")}})
parameter_summary$reliable_mean = sapply(as.numeric(summary$ess_bulk)>100, function(x){if (x){c("Y")} else {c("N")}})
parameter_summary$reliable_quantiles = sapply(as.numeric(summary$ess_tail)>100, function(x){if (x){c("Y")} else {c("N")}})
print(parameter_summary)

write_csv(parameter_summary, "Bayesian_model/testing/data_files/synthetic_parameters.csv")

