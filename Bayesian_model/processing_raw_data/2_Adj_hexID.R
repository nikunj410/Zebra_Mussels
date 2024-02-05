# 1.This code used Adj generated from raw MT data to create a Adj_hex matrix
# 2. Next, we trim Adj_hex to remove hex grids which have less than threshold (here set as 24 = 1 per month) number of incoming and outgoing ships
library(readr)

Adj_raw <- read_csv("Bayesian_model/processing_raw_data/data_files/adj/Adj_raw.csv",show_col_types = FALSE)
Adj_raw = as.matrix(Adj_raw)
rownames(Adj_raw) = colnames(Adj_raw)
portsHexID <- read_csv("Bayesian_model/processing_raw_data/data_files/port_list/portsHexID.csv",show_col_types = FALSE)

N = length(unique(portsHexID$hex_id))

Adj_hex = as.matrix(matrix(0,nrow = N, ncol =N))

rownames(Adj_hex) = 1:N
colnames(Adj_hex) = 1:N

for (i in 1:N) for (j in 1:N) if(i!=j) 
    Adj_hex[i,j] =  sum(Adj_raw[as.character(portsHexID$Port_ID[portsHexID$hex_id == i]),as.character(portsHexID$Port_ID[portsHexID$hex_id == j])])


threshold = 24 # ensuring that each hex grid have at least one incoming and outgoing ship every month 
detele_hex_grd = sort(union(which(colSums(Adj_hex)<threshold),which(rowSums(Adj_hex)<threshold))) # taking the union of incoming and outgoing hex grids that needs to deleted

Adj_hex = Adj_hex[-c(detele_hex_grd),]
Adj_hex = Adj_hex[,-c(detele_hex_grd)]

write_csv(as.data.frame(Adj_hex),"Bayesian_model/processing_raw_data/data_files/adj/Adj_hex.csv")
