library(readr)
library(lamW)       # lambert W function
library(optrees)    # graph operations
library(popbio)
library(reshape2)   # generate edgelist from adjacency matrix
library(rtern)
library(ggplot2)
library(cowplot)
library(doParallel)
registerDoParallel(3)

Leslie_matrix <- function(F2,F3,F4,S1,S2,S3,S4) return(matrix(c(0,S1,0,0,F2,0,S2,0,F3,0,0,S3,F4,0,0,S4), nrow = 4, ncol=4))
de <- function(v,w,B11,B22,B33,B44) return(t(v) %*% diag(c(B11,B22,B33,B44)) %*% w/sum(v*w)) 

arrival_time<- function(F2,F3,F4,S1,S2,S3,S4, B11, B22, B33, B44, gm_Aij){
  Leslie = Leslie_matrix(F2,F3,F4,S1,S2,S3,S4)
  v = reproductive.value(Leslie)
  w = stable.stage(Leslie)
  lambda = lambda(Leslie)
  de = de(v,w,B11, B22,B33,B44)
  return(lambertW0(log(lambda)/(de*gm_Aij))/log(lambda))
} 



port_id = read_csv("Bayesian_model/processing_raw_data/data_files/port_list/portsHexID.csv",show_col_types = FALSE)
Adj = read_csv("Bayesian_model/processing_raw_data/data_files/adj/Adj_hex.csv",show_col_types = FALSE)
Adj <- as.matrix(Adj/sum(Adj))
p = nrow(Adj)
rownames(Adj)=colnames(Adj)
ADJ = as.matrix(Adj)
colnames(ADJ) = 1:p
rownames(ADJ) = 1:p
edgelist = melt(ADJ)
colnames(edgelist) = c("start", "end", "weights")
edgelist = edgelist[which(edgelist$weights != 0),]
edgelist$logTsp = -log(edgelist$weights)
edgelist$weights = -log(edgelist$weights) # Lambert transform distance



fit_zebra = readRDS("Bayesian_model/zebra/data_files/zebra_stan.rds")
draw = as.data.frame(fit_zebra$draws( format="df", variables = c("sigma","f","deff","intro_port")))
chains = max(draw$.chain)
InterPerChain = max(draw$.iteration)
iterations = 1000
draw = draw[,-((ncol(draw)-2):ncol(draw))]

gm_Aij_all = sapply(unique(draw$intro_port),function(x){exp(-mean(spTreeDijkstra(1:p, as.matrix(edgelist[,1:3]), source.node = x, directed = TRUE)$tree.arcs[,3]))})
elasticity_time = as.data.frame(matrix(0, nrow  = iterations, ncol = 12))
colnames(elasticity_time) = c("F2","F3","F4","S1","S2","S3","S4", "B11", "B22", "B33", "B44", "gm_Aij")

OnePlusDelta = 1.001

for(i in 1:iterations)
{
  gm_Aij = gm_Aij_all[which(unique(draw$intro_port) == draw$intro_port[i])]
  F2 = draw$`sigma[1]`[i]*draw$`f[1]`[i]*10^6
  F3 = draw$`sigma[1]`[i]*draw$`f[2]`[i]*10^6
  F4 = draw$`sigma[1]`[i]*draw$`f[3]`[i]*10^6
  S1 = draw$`sigma[2]`[i]
  S2 = draw$`sigma[3]`[i]
  S3 = draw$`sigma[4]`[i]
  S4 = draw$`sigma[5]`[i]
  Leslie = Leslie_matrix(F2,F3,F4,S1,S2,S3,S4)
  v = reproductive.value(Leslie)
  w = stable.stage(Leslie)
  lambda = lambda(Leslie)
  Bii = draw$deff[i]*sum(v*w)/sum(v*w)
  B11 = Bii
  B22 = Bii
  B33 = Bii
  B44 = Bii
  original = arrival_time(F2,F3,F4,S1,S2,S3,S4, B11, B22, B33, B44, gm_Aij)
  
  elasticity_time[i,] = foreach(Npar = 1:12, .combine = 'c') %dopar%{
    -(log(original/arrival_time(Npar == 1 ? (F2*OnePlusDelta) : F2,
                                Npar == 2 ? (F3*OnePlusDelta) : F3,
                                Npar == 3 ? (F4*OnePlusDelta) : F4,
                                Npar == 4 ? (S1*OnePlusDelta) : S1,
                                Npar == 5 ? (S2*OnePlusDelta) : S2,
                                Npar == 6 ? (S3*OnePlusDelta) : S3,
                                Npar == 7 ? (S4*OnePlusDelta) : S4,
                                Npar == 8 ? (B11*OnePlusDelta) : B11,
                                Npar == 9 ? (B22*OnePlusDelta) : B22,
                                Npar == 10 ? (B33*OnePlusDelta) : B33,
                                Npar == 11 ? (B44*OnePlusDelta) : B44,
                                Npar == 12 ? (gm_Aij*OnePlusDelta) : gm_Aij))/log(OnePlusDelta))
  }
}

elasticity_time_ggplot = as.data.frame(cbind(rep(c("F2","F3","F4","S1","S2","S3","S4","B11", "B22", "B33", "B44", "gm_Aij"),iterations),as.vector(t(as.matrix(elasticity_time)))))
colnames(elasticity_time_ggplot) = c("var","elasticity")
elasticity_time_ggplot$elasticity = as.numeric(elasticity_time_ggplot$elasticity)

elasticity_time_ggplot = elasticity_time_ggplot[elasticity_time_ggplot$var != c("gm_Aij"),]
sensitivity  = ggplot(data = elasticity_time_ggplot, aes(x=reorder(var,-elasticity), y=elasticity)) +
  geom_boxplot (colour = "grey50",fill = "grey50",outlier.alpha = 0.01, outlier.size = 0.2)+
  xlab("Parameter")+ylab("Elasticity of Arrival Time")+
  theme( panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(aspect.ratio = 1)+
  scale_y_reverse()+coord_flip()+
  scale_x_discrete(labels = rev(c(expression(S[1]), expression(F[2]),
                                  expression(B[22]), expression(B[11]),
                                  expression(S[2]), expression(F[3]),
                                  expression(B[33]), expression(F[4]),
                                  expression(S[3]), expression(B[44]), expression(S[4]))))

Leslie_draw = as.data.frame(fit_zebra$draws( format="df", variables = "Leslie"))
chains = max(Leslie_draw$.chain)
InterPerChain = max(Leslie_draw$.iteration)
iterations = chains*InterPerChain
Leslie_draw = Leslie_draw[,-((ncol(Leslie_draw)-2):ncol(Leslie_draw))]

set.seed(125)
Nsamples = 100
sample_index = sample(1:iterations, Nsamples, replace=F)
Wmean = colMeans(t(sapply(sample_index, function(x) stable.stage(matrix(unlist(Leslie_draw[x,]), ncol = 4, nrow = 4)))))
Vmean = colMeans(t(sapply(sample_index, function(x) reproductive.value(matrix(unlist(Leslie_draw[x,]), ncol = 4, nrow = 4)))))
relative_dispersal = c(1,1,1,1) ## estimates obtained from Keevin et. al. Journal of Freshwater Ecology (1992)
De_mean = Wmean*relative_dispersal*Vmean
De_mean = De_mean/max(De_mean)

alp = 0.08
scale = 400

WBVgeom_lines <- lapply(sample_index, function(x){
  dat = data.frame(age_group = seq(1:4), WBV = stable.stage(matrix(unlist(Leslie_draw[x,]), ncol = 4, nrow = 4))*
                     relative_dispersal*reproductive.value(matrix(unlist(Leslie_draw[x,]), ncol = 4, nrow = 4)))
  geom_line(data = dat, aes(x= age_group,y = WBV/max(WBV)), alpha = alp ,col="grey50" )})

De = ggplot() + WBVgeom_lines + xlab(c("Age Class")) + ylab("Age-Specific Relative\nContribution to De")+
  scale_x_continuous(breaks = c(1,2,3,4), labels=c("1","2","3","4"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        aspect.ratio = 1)

all_plots = plot_grid(sensitivity, De, ncol = 2, labels = "AUTO", align = "v")
plot(all_plots)

ggsave("Figures/S3_life_history_effects_alternate/sensitivity_alternative.pdf",all_plots, width = 18, height = 8, units = c("cm"))

