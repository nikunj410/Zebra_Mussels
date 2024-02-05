library(lubridate)
library(geodist)
library(readr)


######
Gbif_raw <- read.delim("Bayesian_model/processing_raw_data/data_files/presence_data/gbif.csv")
Gbif_raw = Gbif_raw[-which(Gbif_raw$coordinateUncertaintyInMeters>1000),]

######
NAS_data = read_csv("Bayesian_model/processing_raw_data/data_files/presence_data/NAS.csv",show_col_types = FALSE)
NAS_data = NAS_data[-which(NAS_data$Accuracy == "Approximate"),]
NAS_data$Month[is.na(NAS_data$Month)] = 1
NAS_data$Day[is.na(NAS_data$Day)] = 1
NAS_data= NAS_data[-which(NAS_data$Year>2019),]
NAS_data$eventDate = sapply(1:length(NAS_data$Year), function(x) {paste(as.character(NAS_data$Year[x]),as.character(NAS_data$Month[x]),as.character(NAS_data$Day[x]),sep="-")})
#####

portsHexID <- read_csv("Bayesian_model/processing_raw_data/data_files/port_list/portsHexID.csv",show_col_types = FALSE)
Adj_hex = read_csv("Bayesian_model/processing_raw_data/data_files/adj/Adj_hex.csv",show_col_types = FALSE)
Hex_id = as.integer(colnames(Adj_hex))

presence_dat = data.frame(lat = c(Gbif_raw$decimalLatitude,NAS_data$Latitude),
                          lon = c(Gbif_raw$decimalLongitude,NAS_data$Longitude),
                          date = c(Gbif_raw$eventDate, NAS_data$eventDate),
                          hex_id = NA)

presence_dat$date = decimal_date(as.Date(presence_dat$date))
presence_dat = presence_dat[order(presence_dat$date),]

threshold = 2.5e4 # within 25km radius

distance_matrix = geodist(data.frame(lon = presence_dat$lon, lat = presence_dat$lat),
                          data.frame(lon = portsHexID$Longitude, lat = portsHexID$Latitude),
                          measure = "geodesic")
distance_matrix[which(distance_matrix>threshold)]=NA

presence_dat$hex_id = portsHexID$hex_id[apply(distance_matrix, 1, function(x) which.min(x)[1])]

presence_dat = presence_dat[!is.na(presence_dat$hex_id),]
presence_dat = presence_dat[is.na(match(presence_dat$hex_id,setdiff(presence_dat$hex_id,Hex_id))),]
write_csv(presence_dat,"Bayesian_model/processing_raw_data/data_files/presence_data/zebra_presence.csv")


