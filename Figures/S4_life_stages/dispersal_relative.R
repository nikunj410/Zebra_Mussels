library(popbio)     # calculate w, v, and lambda from Leslie matrix
library(ggplot2)    # plotting
library(cowplot)    # for aligning plots

fit_zebra = readRDS("Bayesian_model/zebra/data_files/zebra_stan.rds")
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
relative_dispersal = c(0,1,1,1) ## estimates obtained from Keevin et. al. Journal of Freshwater Ecology (1992)
De_mean = Wmean*relative_dispersal*Vmean
De_mean = De_mean/max(De_mean)

alp = 0.08
scale = 400
V_col = c("#ED5564")
W_col = c("#4FC1E8")
Wgeom_lines <- lapply(sample_index, function(x){
                        dat = data.frame(age_group = seq(1:4), stable_age = stable.stage(matrix(unlist(Leslie_draw[x,]), ncol = 4, nrow = 4)))
                        geom_line(data = dat, aes(x= age_group,y = stable_age), alpha = alp ,col=W_col )})
Vgeom_lines <- lapply(sample_index, function(x){
                        dat = data.frame(age_group = seq(1:4), rep_value = reproductive.value(matrix(unlist(Leslie_draw[x,]), ncol = 4, nrow = 4)))
                        geom_line(data = dat, aes(x= age_group,y = rep_value/scale), alpha = alp ,col=V_col )})
WBVgeom_lines <- lapply(sample_index, function(x){
                          dat = data.frame(age_group = seq(1:4), WBV = stable.stage(matrix(unlist(Leslie_draw[x,]), ncol = 4, nrow = 4))*
                                             relative_dispersal*reproductive.value(matrix(unlist(Leslie_draw[x,]), ncol = 4, nrow = 4)))
                          geom_line(data = dat, aes(x= age_group,y = WBV/max(WBV)), alpha = alp ,col="grey50" )})

VBW = ggplot() + Wgeom_lines + labs(y = expression(paste("Stable Age Distribution (",bold("w"),")")), color = W_col) + xlab("Age Class")+
  Vgeom_lines + scale_y_continuous(sec.axis = sec_axis(~.*scale, name=expression(paste("Reproductive Value (",bold("v"^"T"),")"))))+
  geom_point(aes(x= seq(1,4), y = Wmean) , col=W_col, size = 3)+
  geom_point(aes(x= seq(1,4), y = Vmean/scale) , col=V_col, size = 3)+
  scale_x_continuous(breaks = c(1,2,3,4), labels=c("1","2","3","4"))+
  theme(axis.text.y  = element_text(colour = W_col), axis.text.y.right  = element_text(colour = V_col),
    axis.title.y = element_text(color = W_col), axis.title.y.right = element_text(color = V_col),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    aspect.ratio = 1) 

plot(VBW)
ggsave("Figures/S4_life_stages/vwcurve.pdf",VBW,   width = 9, height = 7, units = c("cm"))



