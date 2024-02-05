library(readr)      # read files
library(ggplot2)    # plotting
library(cowplot)    # for aligning plots
library(ggrepel) 
library(magrittr)
library(sf)         # plotting USA map
library(ggeasy)

source("Bayesian_model/testing/0_synthetic_parameters.R")
port_id = read_csv("Bayesian_model/processing_raw_data/data_files/port_list/portsHexID.csv",show_col_types = FALSE)
fit_zebra = readRDS("Bayesian_model/zebra/data_files/zebra_stan.rds")
fit_testing = readRDS("Bayesian_model/testing/data_files/synthetic_stan.rds")
draw_zebra = as.data.frame(fit_zebra$draws( format="df", variables = c("intro_year","intro_port")))
draw_testing = as.data.frame(fit_testing$draws( format="df", variables = c("intro_year","intro_port")))
draw_zebra = draw_zebra[,-((ncol(draw_zebra)-2):ncol(draw_zebra))]
draw_testing = draw_testing[,-((ncol(draw_testing)-2):ncol(draw_testing))]

threshold = 1
intro_port_testing <- read_csv("Bayesian_model/testing/data_files/synthetic_intro_port.csv",show_col_types = FALSE)
intro_port_testing = intro_port_testing[intro_port_testing$percent>threshold,]
intro_port_zebra <- read_csv("Bayesian_model/zebra/data_files/intro_port_HexID.csv",show_col_types = FALSE)
intro_port_zebra = intro_port_zebra[intro_port_zebra$percent>threshold,]

Glakes = st_read("Bayesian_model/processing_raw_data/data_files/shape_files/greatlakes",quiet = T) %>% st_transform(crs = 4326)
USA_rivers <- st_read("Bayesian_model/processing_raw_data/data_files/shape_files/NA_rivers",quiet = T) %>% st_transform(crs = 4326)

lakes = data.frame(
  lakes = c( "Lake Erie", "Lake\nSt. Claire"),
  lon = c(-81, -83.1), lat = c( 42.25, 42.6),
  angle = c(25, 0))


p1 = ggplot(data = draw_testing)+
  geom_boxplot(aes(x = intro_year), colour = "grey10",fill = "grey",
               outlier.alpha = 0.1, outlier.size = 1)+
  xlab(expression(tau["10"^-9])) + ylab("") + ylim(-1,1)+ 
  geom_text(aes(label = "*", x = IntroYear, y = 0.1),check_overlap = T, size = 7, fontface = "bold")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())+
  easy_remove_axes(which = c("y"), what = c("line"), teach = FALSE)

p2 = ggplot(data = draw_zebra)+
  geom_boxplot(aes(x = intro_year), colour = "grey10",fill = "grey",
               outlier.alpha = 0.1, outlier.size = 1)+
  xlab(expression(tau["10"^-9])) + ylab("") + ylim(-1,1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())+
  easy_remove_axes(which = c("y"), what = c("line"), teach = FALSE)

p3 = ggplot(intro_port_testing) +
  geom_sf(data = Glakes, color = "white", fill = "white") +
  geom_sf(data = USA_rivers, color = "white") +
  coord_sf(xlim = c(-83.6,-78.8), ylim = c(41.2, 43),expand = FALSE) +
  theme(panel.background = element_rect("grey90"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text(data = port_id, aes(label = "*", x = hex_lon[hex_id==intro_port_hexID][1],
                                y = hex_lat[hex_id==intro_port_hexID][1]),
            check_overlap = T, size = 7, fontface = "bold", nudge_x = 0.25 , nudge_y = -0.2)+
  geom_text(data = lakes, aes(label = lakes, x = lon, y = lat, angle = angle), size = 3.5,
            check_overlap = T, size = 2.5, colour = "grey50")+ xlab("") + ylab("") +
  geom_text_repel(mapping = aes(x = lon, y = lat, label = paste(percent,"%",sep = "")),
                  size = 3, box.padding = 0.5, point.padding = 0, min.segment.length = 0)+
  geom_point( aes(x = lon, y = lat, size = percent), colour = "grey20", alpha = 0.7, show.legend = FALSE)+
  scale_size(range = c(min(intro_port_testing$percent*0.075), max(intro_port_testing$percent*0.075)))+
  scale_y_discrete(breaks = c(43, 42))
p4 = ggplot(intro_port_zebra) +
  geom_sf(data = Glakes, color = "white", fill = "white") +
  geom_sf(data = USA_rivers, color = "white") +
  coord_sf(xlim = c(-83.6,-78.8), ylim = c(41.2, 43),expand = FALSE) +
  theme(panel.background = element_rect("grey90"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text(data = port_id, aes(label = "Kingsville", x = hex_lon[hex_id==intro_port_hexID][1],
                                y = hex_lat[hex_id==intro_port_hexID][1]),
            check_overlap = T, size = 2.5, nudge_x = 0.15, nudge_y = -0.15)+
  geom_text(data = lakes, aes(label = lakes, x = lon, y = lat, angle = angle), size = 3.5, 
            check_overlap = T, size = 2.5, colour = "grey50")+  xlab("") + ylab("") +
  geom_text_repel(mapping = aes(x = lon, y = lat, label = paste(percent,"%",sep = "")),
                  size = 3, box.padding = 0.5, point.padding = 0, min.segment.length = 0)+
  geom_point( aes(x = lon, y = lat, size = percent), colour = "grey20", alpha = 0.7, show.legend = FALSE)+
  scale_size(range = c(min(intro_port_zebra$percent*0.075), max(intro_port_zebra$percent*0.075)))+
  scale_y_discrete(breaks = c(43, 42))

title_left <- ggdraw() + draw_label("Synthetic Data", fontface='bold')
title_right <- ggdraw() + draw_label("Real Data", fontface='bold')
left_col <- plot_grid(title_left, p1, p3, nrow = 3, labels = c("","A","C"), rel_heights = c(0.16,0.75,1.5))
right_col <- plot_grid(title_right, p2, p4, nrow = 3, labels = c("","B","D"), rel_heights = c(0.16,0.75,1.5))

all_plots = plot_grid(left_col, right_col, ncol = 2)
plot(all_plots)

ggsave("Figures/Fig4_intro_port_time/figure3.pdf",all_plots, width = 18 , height = 8, units = c("cm"))
  
