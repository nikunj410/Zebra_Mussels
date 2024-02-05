library(readr)
library(popbio)
library(ggplot2)
library(cowplot)
library(matrixStats)
### posterior predictive checks——first detection time and total observation counts

fit = readRDS("Bayesian_model/testing/data_files/synthetic_stan.rds")
survey = as.matrix(read_csv("Bayesian_model/testing/data_files/synthetic_survey.csv", show_col_types = FALSE))

Survey_draw = as.data.frame(fit$draws( format="df", variables = "Survey"))
chains = max(Survey_draw$.chain)
InterPerChain = max(Survey_draw$.iteration)
iterations = chains*InterPerChain

Survey_draw = Survey_draw[,-((ncol(Survey_draw)-2):ncol(Survey_draw))]
Time = nrow(survey)
n_ports = ncol(survey)
n = ncol(Survey_draw)/(Time*n_ports)
Survey_draw = array(t(Survey_draw),dim = c(n*Time,n_ports,iterations))
dt = 1/3
start_year = 1980

lowerQ = 0.025
upperQ = 0.975

Survey = data.frame(first_detection_obs = sapply(1:n_ports,function(x){which(survey[,x]==1)[1]})*dt+start_year)
Survey$first_detection_predict = rowMedians(sapply(1:iterations,function(y){sapply(1:n_ports,function(x){which(Survey_draw[,x,y]==1)[1]})}), na.rm = T)*dt+start_year
quantile_first_detection= apply(sapply(1:iterations,function(y){sapply(1:n_ports,function(x){which(Survey_draw[,x,y]==1)[1]})})*dt+start_year, 1, quantile, na.rm = T, probs = c(lowerQ,upperQ))
Survey$quantile_first_detection_lower= quantile_first_detection[1,]
Survey$quantile_first_detection_upper= quantile_first_detection[2,]
Survey$first_detection_consistent = as.numeric(Survey$first_detection_obs>Survey$quantile_first_detection_lower & Survey$first_detection_obs<Survey$quantile_first_detection_upper)
Survey$first_detection_consistent = as.factor(Survey$first_detection_consistent)


Survey$mean_counts_obs = colSums(survey)
Survey$mean_counts_predict = rowMedians(sapply(1:iterations, function(x){colSums(Survey_draw[1:Time,,x])}))
quantile_mean_cont = apply(sapply(1:iterations, function(x){colSums(Survey_draw[1:Time,,x])}), 1, quantile, na.rm = T, probs = c(lowerQ,upperQ))
Survey$quantile_mean_cont_lower = quantile_mean_cont[1,]
Survey$quantile_mean_cont_upper = quantile_mean_cont[2,]
Survey$mean_cont_consistent = as.numeric(Survey$mean_counts_obs>Survey$quantile_mean_cont_lower & Survey$mean_counts_obs<Survey$quantile_mean_cont_upper)
Survey$mean_cont_consistent = as.factor(Survey$mean_cont_consistent)


Survey = na.omit(Survey)
Survey = Survey[rev(order(Survey$first_detection_consistent)),]

p1 = ggplot(Survey)+
  geom_function(fun = function(x) x)+ geom_point(aes(x = first_detection_obs, y = first_detection_predict,
                                                     colour = first_detection_consistent)) + 
  scale_color_manual(values = c("0" = "black", "1" = "grey"))+
  xlab("Observed time of\nfirst detection") + ylab ("Predicted time of\nfirst detection")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position = "none")


p2 = ggplot(Survey)+
  geom_function(fun = function(x) x)+ geom_point(aes(x = mean_counts_obs, y = mean_counts_predict,
                                                     colour = mean_cont_consistent) ) + 
  scale_color_manual(values = c("0" = "black", "1" = "grey"))+
  xlab("Number of observed\ndetection events") + ylab ("Number of predicted\ndetection events")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position = "none")
title_left <- ggdraw() + draw_label("Synthetic Data", fontface='bold')
First_detection_syn = plot_grid(title_left,p1, p2, labels = c("","A","C"), align = "v", nrow = 3, rel_heights = c(0.25,3,3))
plot(First_detection_syn)



