# This code clusters the ports into hexagonal grids of 50Km (width). 
# This is done to coarse-grain the shipping network while preserving the large scale structure of the network
library(tidyverse)
library(sf)

width <- 5e4 # clustering ports in hexagonal grids of 50 Km.
ports_dat <- read.csv("Bayesian_model/processing_raw_data/data_files/port_list/portsPlusID.csv")

ports_dat_sf <- st_as_sf(ports_dat, coords = c("Longitude", "Latitude"), crs = 4326) %>% st_transform(2163)
hex_grid <- st_make_grid(ports_dat_sf, c(width, width), square = FALSE) %>% st_as_sf()
hex_grid <- hex_grid[ports_dat_sf, ]
hex_grid <- hex_grid %>% mutate(id = 1:nrow(hex_grid))
ports_dat_sf <- st_join(ports_dat_sf, hex_grid)

hex_centroids <- st_centroid(hex_grid) %>% st_transform(crs = 4326) %>%
  mutate(longitude = st_coordinates(.)[, 1],latitude = st_coordinates(.)[, 2]) %>%
  st_drop_geometry()


ports_dat$hex_id = ports_dat_sf$id
ports_dat$hex_lon = hex_centroids$longitude[match(ports_dat_sf$id,hex_centroids$id)]
ports_dat$hex_lat = hex_centroids$latitude[match(ports_dat_sf$id,hex_centroids$id)]
write_csv(ports_dat,"Bayesian_model/processing_raw_data/data_files/port_list/portsHexID.csv")



