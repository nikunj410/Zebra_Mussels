library(readr)      # read files
library(ggplot2)    # plotting
library(cowplot)    # for aligning plots
library(ggrepel) 

survey = as.matrix(read_csv("Bayesian_model/zebra/data_files/zebra_survey_matrix.csv", show_col_types = FALSE))
n_ports = ncol(survey)
dt = 1/3
start_year = 1980
port_id = read_csv("Bayesian_model/processing_raw_data/data_files/port_list/portsHexID.csv",show_col_types = FALSE)

Adj = read_csv("Bayesian_model/processing_raw_data/data_files/adj/Adj_hex.csv",show_col_types = FALSE)

USA_map <- st_read("Bayesian_model/processing_raw_data/data_files/shape_files/NA_map",quiet = T) %>% st_transform(crs = 4326)
USA_rivers <- st_read("Bayesian_model/processing_raw_data/data_files/shape_files/NA_rivers",quiet = T) %>% st_transform(crs = 4326)
state = c("MO", "IA", "NE", "IL", "MN", "PA", "MI", "KY", "OH", "WV", "IN", "NY", "AL", "MS", "TN", "OK", "AR", "LA", "WI", "KS", "VT")
USA_rivers = USA_rivers[!is.na(match(USA_rivers$state,state)),]
USA_rivers= USA_rivers[is.na(match(USA_rivers$objectid,c(608, 1774, 1772))),]
USA_rivers = USA_rivers[-which(USA_rivers$state=="LA" & USA_rivers$length<50),]


data = data.frame(hex_id = colnames(Adj))
data$colonization_year = as.character(floor(sapply(1:n_ports,function(x){(which(survey[,x]==1)[1]-1)})*dt+start_year))
data$Longitude = port_id$hex_lon[match(colnames(Adj),port_id$hex_id)]
data$Latitude = port_id$hex_lat[match(colnames(Adj),port_id$hex_id)]

intro_year = ggplot(USA_map) + 
  geom_sf(data = USA_map, color = "grey90", fill = "grey90") +
  geom_sf(data = USA_rivers, color = "white") +
  coord_sf(xlim = c(-97,-70.1), ylim = c(28, 49),expand = FALSE) +
  theme(panel.background = element_rect("white"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_point(data = data, aes(x = Longitude, y = Latitude, colour = as.character(colonization_year)))+
  geom_text_repel(data = data, aes(x = Longitude, y = Latitude, label = colonization_year),
                  size = 2, box.padding = 0.25, min.segment.length = 0)+
  guides(colour=guide_legend(ncol=3, title = "Year of First Detection")) +
  theme(legend.position = c(0.82, 0.15),legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA, color = NA),legend.title.align=0.5)

plot(intro_year)
ggsave("Figures/S1_colonization_year/intro_year.pdf",intro_year,   width = 18, height = 18, units = c("cm"))
 
