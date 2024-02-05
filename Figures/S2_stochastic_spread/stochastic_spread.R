library(Matrix)     # for sparce matrix functions
library(readr)      # read files
library(fastmatrix) # vec-permutation operator
library(fBasics)    # vec operator
library(popbio)     # calculate w, v, and lambda from Leslie matrix
library(ggplot2)    # plotting
library(ggpmisc)    # for liner regression in ggplot
library(cowplot)    # for aligning plots
library(magrittr)   # pipes
library(stringr)    # manipulate strings
library(combinat)
library(tictoc)
library(rtern)
library(Hmisc)

Reproduction_Survival = function(pop, A){
  return(c(rpois(1,pop[2]*A[1,2]+pop[3]*A[1,3]+pop[4]*A[1,4]),
           rbinom(2,pop[1:2],c(A[2, 1],A[3, 2])),
           sum(rbinom(2,pop[3:4],A[4,3:4]))))}


lowerQ = 0.025
upperQ = 0.975

arrival_theory = read_csv("Figures/Fig1_metapopulation/arrival_stats_theory.csv",show_col_types = FALSE)
arrival_theory = arrival_theory[,c(1,2,14)]
port_id = read_csv("Bayesian_model/processing_raw_data/data_files/port_list/portsHexID.csv",show_col_types = FALSE)
Adj = read_csv("Bayesian_model/processing_raw_data/data_files/adj/Adj_hex.csv",show_col_types = FALSE)
Adj <- as.matrix(Adj/sum(Adj))
rownames(Adj)=colnames(Adj)
p = nrow(Adj)
presence_data = read_csv("Bayesian_model/processing_raw_data/data_files/presence_data/zebra_presence.csv",show_col_types = FALSE)

million = 1e6
s = 4
Leslie = Matrix(matrix(rep(0,s^2), nrow = s), sparse =T)
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

hr = 1/3
w = stable.stage(Leslie)
v = reproductive.value(Leslie)
lambda = lambda(Leslie)
disper_ages = c(1,1,1,1)
De = hr*sum(v*disper_ages*w)/sum(v*w)

source = which(rownames(Adj)==presence_data$hex_id[1])
Intro_port_name = tail(port_id$Name[which(port_id$hex_id==presence_data$hex_id[1])],1)

#Permutation (P) (Hunter and Caswell, Ecolological Modelling (2005))
P = commutation(m = s, n = p)
P = sparseMatrix(P$row,P$col,x=1)

redistribution_matrix = array(rep(0,p*p*s), dim = c(p,p,s))
for (i in 1:s) redistribution_matrix[,,i] = t(diag(p)-(diag(colSums(t(Adj)))-t(Adj))*disper_ages[i]*hr)

K = rep(ceiling(1e9*w), p) # setting carrying capacity because R cannot handle integers above 1e10
K = P%*%K
set.seed(123)
## this is where loop happens
tic()
iteration = 0
while(iteration < 1000){
  print(iteration)
  popN = Matrix(matrix(rep(0,s*p),nrow = s, ncol =p), sparse=T)
  initial_pop = rep(0,s)

  while(sum(initial_pop)==0) initial_pop = rpois(4,w)
  popN[,source] = initial_pop
  pop_stages = Matrix(fBasics::vec(t(popN)), sparse=T)
  popW = (matrix(pop_stages, nrow = p, ncol = s)%*%v)/sum(v*w)

  arrival_stats = as.data.frame(matrix(rep(0,2*p),nrow=p,ncol=2))
  colnames(arrival_stats)=c("port_Hex_id","time_arrival")
  invaded = 1
  arrival_stats[1,]=c(source,0)
  time=0
  flag = 0
  while (invaded<p){
    time=time+1

    pop_stages_disp = pop_stages
    pop_stages = pmin(pop_stages,K) # setting carrying capacity because R cannot handle integers above 1e10
    # pop_stages_disp = as.vector(sapply(1:s, function(i){rowSums(sapply(1:p,function(j){rmultinom(1,pop_stages[(i-1)*p+j], redistribution_matrix[j, ,i])}))}))
    pop_stages_disp = as.vector(sapply(1:s, function(i){sum(pop_stages[((1-p):0)+p*i])>0 ? rowSums(sapply(which(pop_stages[((1-p):0)+p*i]>0),function(j){rmultinom(1,pop_stages[(i-1)*p+j], redistribution_matrix[j,,i])})):rep(0,p)}))
    
    pop_stages_disp = t(P)%*%pop_stages_disp
    pop_stages_rep_sur = as.vector(sapply(1:p,function(i){sum(pop_stages_disp[((i-1)*s+1):(s*i)])>0 ? Reproduction_Survival(pop_stages_disp[((i-1)*s+1):(s*i)], Leslie) : rep(0,s)}))
    pop_stages = P%*%pop_stages_rep_sur

    popW = (matrix(pop_stages, nrow = p, ncol = s)%*%v)/sum(v*w)
    if (sum(popW)==0) { flag = 1; break}
    invaded_ports = which(popW>=1)
    new_invasion = setdiff(invaded_ports, arrival_stats$port_Hex_id)
    if (length(new_invasion)!=0){
      arrival_stats[(invaded+1):(invaded+length(new_invasion)),] = cbind(new_invasion[],time-log(popW[new_invasion])/log(lambda))
      invaded=invaded+length(new_invasion)
    }
  }
  
  if(flag == 0){
    arrival_stats = arrival_stats[order(arrival_stats$time_arrival),]
    arrival_stats$port_Hex_id = as.numeric(colnames(Adj)[arrival_stats$port_Hex_id])
    arrival_stats = arrival_stats[match(arrival_theory$port_Hex_id,arrival_stats$port_Hex_id),]
    time_stochastic = arrival_stats$time_arrival
    arrival_theory = cbind(arrival_theory, time_stochastic)
    colnames(arrival_theory)[ncol(arrival_theory)] = c(paste("iteration_",iteration+1,sep=""))
    iteration = iteration+1
  }
  else {flag == 0}
}
toc()
write.csv(arrival_theory, "Figures/S2_stochastic_spread/arrival_stats_stochastic.csv", row.names=FALSE)

comparision = data.frame(port_Hex_id = arrival_theory$port_Hex_id,
                         theory = arrival_theory$true_time_approx2,
                         deterministic = arrival_theory$time_arrival,
                         stochastic = rowMeans(as.matrix(arrival_theory[,4:ncol(arrival_theory)])))
quantile = apply(as.matrix(arrival_theory[,4:ncol(arrival_theory)]), 1, quantile, na.rm = T, probs = c(lowerQ,upperQ))
comparision$quantile_lower = quantile[1,]
comparision$quantile_upper = quantile[2,]

p1 = ggplot(comparision, aes(y = theory, x = deterministic,
                             ymin = deterministic,
                             ymax = deterministic))+
  geom_function(fun = function(x) x)+ geom_pointrange(fill='black', color='grey', shape=21, fatten = 3, alpha = 0.8) + 
  ylab("Theoretical\nExpectation") + xlab ("Deterministic\nSimulation")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p2 = ggplot(comparision, aes(x = deterministic, y = stochastic,
                        ymin = quantile_lower,
                        ymax = quantile_upper))+
  geom_function(fun = function(x) x)+ geom_pointrange(fill='black', color='grey', shape=21, fatten = 3, alpha = 0.8) + 
  xlab("Deterministic\nSimulation") + ylab ("Stochastic\nSimulation")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

comp = plot_grid(p1, p2, labels = c("A","B"), ncol = 2, align = "h", axis = "b", rel_heights  = c(1,1,1))
plot(comp)
ggsave("Figures/S2_stochastic_spread/comparison.pdf",comp, width = 14, height = 7, units = c("cm"))


