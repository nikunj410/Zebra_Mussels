library(sp)
portsHexID <- read_csv("Bayesian_model/processing_raw_data/data_files/port_list/portsHexID.csv",show_col_types = FALSE)
hull = list(lat = c(43.1806,42.2158,40.9031,42.1716), lon = c(-78.5357,-78.4380,-83.2566,-84.1128) )
Adj = read_csv("Bayesian_model/processing_raw_data/data_files/adj/Adj_hex.csv",show_col_types = FALSE)

prior_intro_port = unique(portsHexID$hex_id[as.logical(point.in.polygon(portsHexID$Longitude, portsHexID$Latitude,
                                                     hull$lon, hull$lat, mode.checked=FALSE))])
prior_intro_port = intersect(prior_intro_port,as.integer(colnames(Adj)))
write_csv(as.data.frame(prior_intro_port),"Bayesian_model/processing_raw_data/data_files/intro_port_prior/intro_port_prior_hexID.csv")
