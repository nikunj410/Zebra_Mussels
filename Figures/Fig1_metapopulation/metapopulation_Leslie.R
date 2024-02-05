library(Matrix)     # for sparce matrix functions
library(readr)      # read files
library(fastmatrix) # vec-permutation operator
library(fBasics)    # vec operator
library(lamW)       # lambert W function
library(optrees)    # graph operations
library(popbio)     # calculate w, v, and lambda from Leslie matrix
library(reshape2)   # generate edgelist from adjacency matrix
library(sf)         # plotting USA map
library(ggforce)    # draw circle
library(ggplot2)    # plotting
library(ggpmisc)    # for liner regression in ggplot
library(cowplot)    # for aligning plots
library(magrittr)   # pipes
library(stringr)    # manipulate strings
library(igraph)

gc()
USA_map <- st_read("Bayesian_model/processing_raw_data/data_files/shape_files/NA_map",quiet = T) %>% st_transform(crs = 4326)
USA_rivers <- st_read("Bayesian_model/processing_raw_data/data_files/shape_files/NA_rivers",quiet = T) %>% st_transform(crs = 4326)
state = c("MO", "IA", "NE", "IL", "MN", "PA", "MI", "KY", "OH", "WV", "IN", "NY", "AL", "MS", "TN", "OK", "AR", "LA", "WI", "KS", "VT")
USA_rivers = USA_rivers[!is.na(match(USA_rivers$state,state)),]
USA_rivers= USA_rivers[is.na(match(USA_rivers$objectid,c(608, 1774, 1772))),]
USA_rivers = USA_rivers[-which(USA_rivers$state=="LA" & USA_rivers$length<50),]

port_id = read_csv("Bayesian_model/processing_raw_data/data_files/port_list/portsHexID.csv",show_col_types = FALSE)
Adj = read_csv("Bayesian_model/processing_raw_data/data_files/adj/Adj_hex.csv",show_col_types = FALSE)
Adj <- as.matrix(Adj/sum(Adj))
rownames(Adj)=colnames(Adj)
Adj = Matrix(Adj, sparse=T)
p = nrow(Adj)
presence_data = read_csv("Bayesian_model/processing_raw_data/data_files/presence_data/zebra_presence.csv",show_col_types = FALSE)

river_distance = list(Hex_id = as.integer(colnames(Adj)),
                      distance = sapply(1:ncol(Adj), function(x) min(port_id$RivDist_Kingsville_km[which(as.integer(colnames(Adj))[x]==port_id$hex_id)])))

million = 1e6
s = 4
Leslie = Matrix(matrix(rep(0,s^2), nrow = s), sparse =T)
# Parameters of the Leslie Matrix were taken from Casagrandi et. al. Freshwater Biology (2007)
sigma = c(0.014, 0.61, 0.49, 0.25, 0.087)
f = c(0.14, 0.23, 0.5)
Leslie[1, 2] = sigma[1] * f[1] * million
Leslie[1, 3] = sigma[1] * f[2] * million
Leslie[1, 4] = sigma[1] * f[3] * million
Leslie[2, 1] = sigma[2]
Leslie[3, 2] = sigma[3]
Leslie[4, 3] = sigma[4]
Leslie[4, 4] = sigma[5]

dt = 1/3 # time resolution is four months = 1/3 year
Leslie = Leslie^dt # re-scaling Leslie matrix for a 4 month cycle

hr = 1/3
w = stable.stage(Leslie)
v = reproductive.value(Leslie)
lambda = lambda(Leslie)
disper_ages = c(1,1,1,1)
De = hr*sum(v*disper_ages*w)/sum(v*w)

source = which(rownames(Adj)==presence_data$hex_id[1])

Intro_port_name = tail(port_id$Name[which(port_id$hex_id==presence_data$hex_id[1])],1)

#Permutation (P), block demography (B) and block dispersal (M) Matrix (Hunter and Caswell, Ecolological Modelling (2005)) 
P = commutation(m = s, n = p)
P = sparseMatrix(P$row,P$col,x=1)


B = bdiag(replicate(p,Leslie,simplify=FALSE))
M = bdiag(sapply(1:s, function(x){diag(p)-(diag(colSums(t(Adj)))-t(Adj))*disper_ages[x]*hr}))

# Dispersal first then demography (Hunter and Caswell, Ecolological Modelling (2005)) 
A = P%*%B%*%t(P)%*%M

popN = Matrix(matrix(rep(0,s*p),nrow = s, ncol =p), sparse=T)
popN[,source] = w

pop_stages = Matrix(fBasics::vec(t(popN)), sparse=T)
popW = (matrix(pop_stages, nrow = p, ncol = s)%*%v)/sum(v*w)

arrival_stats = as.data.frame(matrix(rep(0,2*p),nrow=p,ncol=2))
colnames(arrival_stats)=c("port_Hex_id","time_arrival")
invaded = 1
arrival_stats[1,]=c(source,0)
time=0

while (invaded<p){
  time=time+1
  pop_stages = A%*%pop_stages
  popW = (matrix(pop_stages, nrow = p, ncol = s)%*%v)/sum(v*w)
  invaded_ports = which(popW>=1)
  new_invasion = setdiff(invaded_ports, arrival_stats$port_Hex_id)
  if (length(new_invasion)!=0){
    arrival_stats[(invaded+1):(invaded+length(new_invasion)),] = cbind(new_invasion[],time-log(popW[new_invasion])/log(lambda))
    invaded=invaded+length(new_invasion)
  }
}

arrival_stats = arrival_stats[order(arrival_stats$time_arrival),]

ADJ = as.matrix(Adj)
colnames(ADJ) = 1:p
rownames(ADJ) = 1:p
edgelist = melt(ADJ)
colnames(edgelist) = c("start", "end", "weights")
edgelist = edgelist[which(edgelist$weights != 0),]
edgelist$logTsp = -log(edgelist$weights)
edgelist$weights = lambertW0(log(lambda)/(De*edgelist$weights)) # Lambert transform distance
spTree = spTreeDijkstra(1:p, as.matrix(edgelist[,1:3]), source.node = source, directed = TRUE)
spTree_logTsp = spTreeDijkstra(1:p, as.matrix(edgelist[,c(1,2,4)]), source.node = source, directed = TRUE)

NetDistM = ADJ
NetDistM[which(NetDistM!=0)]= -log(NetDistM[which(NetDistM!=0)])
NetDistM_graph = graph_from_adjacency_matrix(as.matrix(NetDistM),weighted=T)

invasion_path = shortest_paths(NetDistM_graph, from = source, to = V(NetDistM_graph),
                               mode = c("out"), output = c("vpath"),
                               algorithm = c("dijkstra"))



arrival_stats$river_distance = river_distance$distance[spTree$tree.nodes]
arrival_stats$network_dist = spTree$distances[spTree$tree.nodes]
arrival_stats$flow = 0.5*(colSums(Adj)+rowSums(Adj))[spTree$tree.nodes]*1.55
arrival_stats$geo_lat = port_id$hex_lat[match(colnames(Adj), port_id$hex_id)][spTree$tree.nodes]
arrival_stats$geo_lon = port_id$hex_lon[match(colnames(Adj), port_id$hex_id)][spTree$tree.nodes]
arrival_stats$port_Hex_id = colnames(Adj)[spTree$tree.nodes]


shortest_tree = as.data.frame(spTree$tree.arcs)
shortest_graph = graph(as.vector(t(as.matrix(shortest_tree[,1:2]))))
shortest.layout = layout_as_tree(shortest_graph, root = source, circular=T)
shortest.layout[,2] = shortest.layout[,2]-shortest.layout[source,2]
shortest.layout = spTree$distances*(shortest.layout/sqrt(rowSums(shortest.layout^2)))
shortest.layout[source,]=c(0,0)

arrival_stats$port_Hex_id = colnames(Adj)[spTree$tree.nodes]
arrival_stats$parent_Hex_id = colnames(Adj)[c(source,spTree$tree.arcs[,1])]
arrival_stats$net_lat = shortest.layout[spTree$tree.nodes,2]
arrival_stats$net_lon = shortest.layout[spTree$tree.nodes,1]


circles <- data.frame(
  x0 = rep(0, 5), y0 = rep(0, 5), 
  r = seq(0, 10*max(arrival_stats$network_dist)%/%10,
          length.out = max(arrival_stats$network_dist)%/%10+1))


p1 = ggplot(arrival_stats) + 
  geom_smooth(aes(x = river_distance, y = time_arrival),method = "lm",
              formula=y~0+x, se = F, colour = "grey", size = 1)+
  geom_point(aes(x = river_distance, y = time_arrival, colour = time_arrival, size = flow),
             alpha = 0.5, show.legend = FALSE) +
  scale_size(range = c(0.5,8))+
  stat_poly_eq(aes(x = river_distance, y = time_arrival),formula=y~0+x)+
  geom_point(aes(x = 0, y = 0), size = 2, colour = "saddlebrown")+
  geom_text(aes(label = "Kingsville" , x = 600, y = 0),check_overlap = T, size = 4, colour = "saddlebrown" )+
  xlab("River Distance (Km)") + ylab("Simulated Arrival Time") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio=1)+
  theme(legend.background = element_rect(fill = NA),legend.key = element_rect(fill = NA, color = NA),
        legend.position = c(0.85, 0.3),legend.key.size = unit(0.6, "cm"),legend.key.width = unit(0.5,"cm"),
        legend.title.align=0.5)

p2 = ggplot(USA_map) + 
  geom_sf(data = USA_map, color = "grey90", fill = "grey90") +
  geom_sf(data = USA_rivers, color = "white") + 
  coord_sf(xlim = c(-97,-70.1), ylim = c(28, 49),expand = FALSE) +
  theme(panel.background = element_rect("white"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_segment(data = arrival_stats, aes(xend = geo_lon, yend = geo_lat,
                                         x = geo_lon[match(parent_Hex_id,port_Hex_id)],
                                         y = geo_lat[match(parent_Hex_id,port_Hex_id)]),
               arrow = arrow(length = unit(0.075,"cm"), type = "closed"), colour = "black", size = 0.1)+
  geom_point(data = arrival_stats, aes(x = geo_lon, y = geo_lat,colour = time_arrival,
                                       size = flow),show.legend = FALSE, alpha = 0.5)+
  scale_size(range = c(0.5,8))+
  geom_point(data = arrival_stats, aes(x = head(geo_lon,1), y = head(geo_lat,1)),
             size = 2, colour = "saddlebrown")+ xlab("") + ylab("")


p3 = ggplot(arrival_stats) + 
  geom_smooth(aes(x = network_dist,y = time_arrival), method = "lm",
              formula=y~0+x, se = F, colour = "grey", size = 1)+
  stat_poly_eq(aes(x = network_dist, y = time_arrival),formula=y~0+x)+
  geom_point(aes(x = network_dist,y = time_arrival, colour = time_arrival,
                 size = flow), alpha = 0.5) +
  guides(size =  guide_legend(title = "Shipping\nFlux"),
         colour = guide_legend(title = "Arrival\nTime", override.aes = list(alpha = 0.5, size = 4)))+
  scale_size_continuous(breaks = c(.01, .04, 0.16), labels = c("1k", "4k", "16k"), range = c(0.5, 8))+
  geom_point(aes(x = 0, y = 0), size = 2, colour = "saddlebrown")+
  xlab("Network Distance") + ylab("") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio=1)+
  theme(legend.background = element_rect(fill = NA),legend.key = element_rect(fill = NA, color = NA),
        legend.position = c(0.8, 0.3),legend.key.size = unit(0.6, "cm"),legend.key.width = unit(0.2,"cm"),
        legend.title.align=0.5, legend.box = "horizontal")

p4 = ggplot(arrival_stats, aes(x = net_lon,y = net_lat)) + 
  geom_circle(aes(x0 = x0, y0 = y0, r = r), data = circles, size = 1, colour = "white", inherit.aes = FALSE)+
  geom_segment(aes(xend = net_lon, yend = net_lat, x = net_lon[match(parent_Hex_id,port_Hex_id)],
                   y = net_lat[match(parent_Hex_id,port_Hex_id)]), colour = "black", size = 0.1,
               arrow = arrow(length = unit(0.075,"cm"), type = "closed"))+
  geom_point(aes(colour = time_arrival,size = flow),show.legend = FALSE, alpha = 0.5)+
  scale_size(range = c(0.5,8))+
  geom_point(aes(x = 0, y = 0), size = 2, colour = "saddlebrown") + xlab("") + ylab("") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), aspect.ratio=1)

title_net <- ggdraw() + draw_label("Network Space", fontface='bold')
title_geo <- ggdraw() + draw_label("Geographical Space", fontface='bold')
geo_plots = plot_grid(title_geo, p1, p2, nrow = 3, labels = c("","A","C"), rel_heights = c(0.75,10,10))
net_plots = plot_grid(title_net, p3, p4, nrow = 3, labels = c("","B","D"), rel_heights = c(0.75,10,10))
all_plots = plot_grid(geo_plots, net_plots, ncol = 2,  align = "v")

plot(all_plots)
ggsave("Figures/Fig1_metapopulation/figure1.pdf",all_plots,   width = 18, height = 18, units = c("cm"))


arrival_stats$true_time_approx1 = spTree$distances[spTree$tree.nodes]/log(lambda)
for (i in 1:p){
  nsp = 0
  child_id = arrival_stats$port_Hex_id[i] 
  parent_id = c("999")
  while(parent_id != c("129") & child_id != c("129")){
    parent_id = arrival_stats$parent_Hex_id[which(arrival_stats$port_Hex_id == child_id)]
    nsp = nsp+1
    child_id = parent_id}
  arrival_stats$nsp[i] = nsp 
}
arrival_stats$Tsp = exp(-spTree_logTsp$distances[spTree_logTsp$tree.nodes])
 

### deviation between shortest network path and most probable network path
NetDistM = ADJ
NetDistM[which(NetDistM!=0)]= -log(NetDistM[which(NetDistM!=0)])
NetDistM_graph = graph_from_adjacency_matrix(as.matrix(NetDistM),weighted=T)

invasion_path = shortest_paths(NetDistM_graph, from = source, to = V(NetDistM_graph),
                               mode = c("out"), output = c("vpath"),
                               algorithm = c("dijkstra"))

nsp = sapply(1:p, function(i) length(invasion_path$vpath[[i]]))-1
Prod_Pij = exp(-distances(NetDistM_graph, v = source, to = V(NetDistM_graph), mode = c("out"), algorithm = c( "dijkstra")))
arrival_stats$true_time_approx2 = c(0,sort(sapply((1:p)[-c(source)],
                                   function(i) (nsp[i]/log(lambda))*lambertW0((log(lambda)/(nsp[i]*De))*(factorial(nsp[i])/Prod_Pij[i])^(1/nsp[i])))))
hist(100*abs(arrival_stats$time_arrival[2:p]-arrival_stats$true_time_approx2[2:p])/arrival_stats$time_arrival[2:p], xlab = "% error", main = "")
percent_error = mean(abs(arrival_stats$time_arrival[2:p]-arrival_stats$true_time_approx2[2:p])/arrival_stats$time_arrival[2:p])*100

write.csv(arrival_stats, "Figures/Fig1_metapopulation/arrival_stats_theory.csv", row.names=FALSE)
