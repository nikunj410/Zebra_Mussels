library(geodist)
library(readr)
portsHexID <- read_csv("Bayesian_model/processing_raw_data/data_files/port_list/portsHexID.csv",show_col_types = FALSE)
Adj_hex = read_csv("Bayesian_model/processing_raw_data/data_files/adj/Adj_hex.csv",show_col_types = FALSE)
Hex_id = as.integer(colnames(Adj_hex))

Hex_grd = data.frame(hex_id = unique(portsHexID$hex_id), lat = unique(portsHexID$hex_lat), lon = unique(portsHexID$hex_lon))
Hex_grd = Hex_grd[order(Hex_grd$hex_id),]
Hex_grd = Hex_grd[-setdiff(Hex_grd$hex_id,Hex_id),]

DistM = geodist(Hex_grd, measure = "geodesic")
DistM = as.data.frame(DistM/sd(DistM[DistM !=0]))

rownames(DistM) = Hex_id
colnames(DistM) = Hex_id

write_csv(DistM,"Bayesian_model/processing_raw_data/data_files/DistM/DistM.csv")


