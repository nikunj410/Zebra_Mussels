library(popbio)
library(Matrix)
library(boot)
library(MASS)

Adj = read_csv("Bayesian_model/processing_raw_data/data_files/adj/Adj_hex.csv",show_col_types = FALSE)
DistM = read_csv("Bayesian_model/processing_raw_data/data_files/DistM/DistM.csv",show_col_types = FALSE)
N = ncol(Adj)

Time = 120
intro_port_hexID = 129
s = 4
K = 1e12
No = 1e4;

p = 5e-13
Eta = 1.5 
Rho = 0.5
Sigma = 0.75
beta_mu = -2

H = K*p
P_Kd = 1-exp(-H)
P_mus = inv.logit(beta_mu)
N_knot = No/K
start_year = 1980
N_star = 1e-10


Leslie = Matrix(matrix(rep(0,s^2), nrow = s), sparse =T)
million = 1e6
# Parameters of the Leslie Matrix were taken from Casagrandi et. al. Freshwater Biology (2007)
sigma = c(0.014, 0.61, 0.49, 0.25, 0.087)
f = c(0.14, 0.23, 0.5)
Leslie[1, 2] = sigma[1] * f[1] * million
Leslie[1, 3] = sigma[1] * f[2] * million
Leslie[1, 4] = sigma[1] * f[3] * million
Leslie[2, 1] = sigma[2]
Leslie[3, 2] = sigma[3]
Leslie[4, 3] = sigma[4]
Leslie[4, 4] = sigma[5]

dt = 1/3 # time resolution is four months = 1/3 year
Leslie = Leslie^dt # re-scaling Leslie matrix for a 4 month cycle

deff = 1/3
w = stable.stage(Leslie)
v = reproductive.value(Leslie)
Lambda = lambda(Leslie)
r = log(Lambda)/dt
IntroYear = log(N_star/N_knot)/r+start_year
Deff = deff/dt
disper_ages = c(1,1,1,1)
hr = deff/sum(v*disper_ages*w)/sum(v*w)

set.seed(123)

SIGMA = (Eta^2)*exp(-(DistM^2)/(2*Rho^2)) + (sigma^2)*diag(N)
p_sampling  = inv.logit(as.vector(mvrnorm(n = 1, mu = rep(beta_mu,N), Sigma = SIGMA)))
