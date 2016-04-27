library(igraph)
library(stringdist)
library(dplyr)

### 0. ent resolution ----
#similarity index
simmilarity_index <- function(x, y){
  
  sim_index <- 1- stringdist(x, y, method = "lv")/max(nchar(as.character(x)), nchar(as.character(y)))
  return(sim_index)
  
}

#entity resolution
edge_list <- mutate_each(edge_list, funs(tolower))
edge_list <- as.data.frame(sapply(edge_list, function(x) gsub("latin1", "", x)))
edge_list <- as.data.frame(sapply(edge_list, function(x) gsub("ASCII", "", x)))
edge_list <- as.data.frame(sapply(edge_list, function(x) gsub('@stategov|@stategoy', "", x)))
#edge_list <- as.data.frame(sapply(edge_list, function(x) gsub("[[:punct:]]", "", x)))
edge_list <- as.data.frame(sapply(edge_list, function(x) gsub('hdr@clintonemailcom', "h", x)))


for (i in 1:nrow(node_names)){
  for (j in 1:i){
    if (simmilarity_index(edge_list$from[i],edge_list$from[j]) > 0.7){
      edge_list$from[i] <- edge_list$from[j]
    }
  }
}

for (i in 1:nrow(node_names)){
  for (j in 1:i){
    if (simmilarity_index(edge_list$V1[i],edge_list$V1[j]) > 0.7){
      edge_list$V1[i] <- edge_list$V1[j]
    }
  }
}

write.csv(edge_list, file = "edge_list.csv")

### 1. importing data ----
#  Import  edge list with classifications
edge_classification <- read.csv('edge_list.csv', sep = ",")
edge_classification$X <- NULL

#seperating clsasified and non-classified edges
classified <- edge_classification[ is.na(edge_classification$reason),]
classified <- classified[c(1,2)]

unclassified <- edge_classification[ !is.na(edge_classification$reason),]
unclassified <- unclassified[c(1,2)]

### 2. create adjacency matrix ----
graphized <- function(x){
  mat <- as.matrix(get.adjacency(graph.edgelist(as.matrix(x), directed=T)))
  graph <- graph.adjacency(mat,mode="directed", weighted = TRUE)
  return(graph)
}

graph_classified <- graphized(classified)
graph_unclassified <- graphized(unclassified)
#node_names <- as.data.frame(V(hillary_graph)$name)

### 3. overall topography of the network ----
topo <- function(x){
  nodes <- vcount(x)
  edges <- ecount(x)
  density <- graph.density(x, loops = T)
  
  topo <- data.frame(nodes,edges,density)
  return(topo)
}

topo_classified <- topo(graph_classified)
topo_unclassified <- topo(graph_unclassified)

### 4. centrality measures ----

cntrlty_measures <- function(x){
  inDegreeC <- degree(x, mode="in", loops = TRUE, normalized = FALSE)
  outDegreeC <- degree(x, mode="out", loops = TRUE, normalized = FALSE)
  totalDegreeC <- degree(x)
  inClosenessC <- closeness(graph_classified, mode='in')
  outClosenessC <- closeness(graph_classified, mode='out')
  totalClosenessC <- closeness(graph_classified)
  betweennessC <- betweenness(graph_classified, directed = T)
  eigenC <- evcent(graph_classified)
  #bonC <- bonpow(graph_classified)
  
  s <- data.frame(inDegreeC, outDegreeC, totalDegreeC, inClosenessC, outClosenessC, totalClosenessC, betweennessC, eigenC)
  s = s[,c(1:8)]
  
  return(s)
}

centrality_classified <- cntrlty_measures(graph_classified)

#for unclassified network
inDegreeU <- degree(graph_unclassified, mode="in", loops = TRUE, normalized = FALSE)
outDegreeU <- degree(graph_unclassified, mode="out", loops = TRUE, normalized = FALSE)
totalDegreeU <- degree(graph_unclassified)
inClosenessU <- closeness(graph_unclassified, mode='in')
outClosenessU <- closeness(graph_unclassified, mode='out')
totalClosenessU <- closeness(graph_unclassified)
betweennessU <- betweenness(graph_unclassified, directed = T)
eigenU <- evcent(graph_unclassified)
#bonU <- bonpow(graph_unclassified)

centrality_unclassified <- data.frame(inDegreeU, outDegreeU, totalDegreeU, inClosenessU, outClosenessU, totalClosenessU, betweennessU, eigenU)
centrality_unclassified = centrality_unclassified[,c(1:8)]

# correlate the measures
corClassified <- cor(centrality_classified)
corUn <- cor(centrality_unclassified)

