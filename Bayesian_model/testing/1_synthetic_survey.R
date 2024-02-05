library(cmdstanr)
library(readr)
library(MASS)
library(Rlab)
library(boot)

options(mc.cores = parallel::detectCores())
source("Bayesian_model/testing/0_synthetic_parameters.R")
Adj = read_csv("Bayesian_model/processing_raw_data/data_files/adj/Adj_hex.csv",show_col_types = FALSE)
DistM = read_csv("Bayesian_model/processing_raw_data/data_files/DistM/DistM.csv",show_col_types = FALSE)

A <- as.matrix(Adj/sum(Adj))
rownames(A)=colnames(A)
N = length(rownames(A))

intro_index = which(colnames(A)==intro_port_hexID) 
survey_sim = matrix(0,nrow = Time, ncol = N)


set.seed(123)

SIGMA = (Eta^2)*exp(-(DistM^2)/(2*Rho^2)) + (sigma^2)*diag(N)
p_sampling  = inv.logit(as.vector(mvrnorm(n = 1, mu = rep(beta_mu,N), Sigma = SIGMA)))
Q_tran = t(A - diag(rowSums(A)));
P_tran = deff*Q_tran+diag(N)
Pop = matrix(0,nrow = N, ncol = Time)

for (i in 1:Time){
  if(i==1) {
    Pop[,i] = array(0,N)
    Pop[intro_index] = No
  }
  else{
    growth = Lambda^(1-Pop[,i-1]/K)
    Pop[,i] = P_tran%*%(growth*Pop[,i-1])
  }
  survey_sim[i,] = rbern(N,p_sampling*(1-exp(-Pop[,i]*p)))
}

Pop = Pop/K
colnames(Pop) = c(1:Time)
rownames(Pop) = colnames(Adj)
colnames(survey_sim) = colnames(Adj)
write_csv(as.data.frame(survey_sim),"Bayesian_model/testing/data_files/synthetic_survey.csv")
print(sum(survey_sim))
