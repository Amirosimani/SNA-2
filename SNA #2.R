library(igraph)
library(stringdist)
library(dplyr)

### ent resolution ----
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

### importing data ----
#  Import  edge list with classifications
edge_classification <- read.csv('edge_list.csv', sep = ",")
edge_classification$X <- NULL

#seperating clsasified and non-classified edges
classified <- edge_classification[ is.na(edge_classification$reason),]
classified <- classified[c(1,2)]

unclassified <- edge_classification[ !is.na(edge_classification$reason),]
unclassified <- unclassified[c(1,2)]

### creat adjacency matrix ----
mat_classified <- as.matrix(get.adjacency(graph.edgelist(as.matrix(classified), directed=T)))
mat_unclassified <- as.matrix(get.adjacency(graph.edgelist(as.matrix(unclassified), directed=T)))

# Create igraph object from this matrix
graph_classified <- graph.adjacency(mat_classified,mode="directed", weighted = TRUE)
graph_unclassified <- graph.adjacency(mat_unclassified,mode="directed", weighted = TRUE)

#node_names <- as.data.frame(V(hillary_graph)$name)

### overall topography of the network ----
vcount(hillary_graph)
ecount(hillary_graph)

graph.density(hillary_graph, loops=T)

### centrality measures ----
#for classified network "turn it to a function"
inDegreeC <- degree(graph_classified, mode="in", loops = TRUE, normalized = FALSE)
outDegreeC <- degree(graph_classified, mode="out", loops = TRUE, normalized = FALSE)
totalDegreeC <- degree(graph_classified)
inClosenessC <- closeness(graph_classified, mode='in')
outClosenessC <- closeness(graph_classified, mode='out')
totalClosenessC <- closeness(graph_classified)
betweennessC <- betweenness(graph_classified, directed = T)
eigenC <- evcent(graph_classified)
bonC <- bonpow(graph_classified)

sumC <- data.frame(inDegreeC, outDegreeC, totalDegreeC, inClosenessC, outClosenessC, totalClosenessC, betweennessC, bonC, eigenC)
sumC = sumC[,c(1:9)]

#for unclassified network
inDegreeU <- degree(graph_unclassified, mode="in", loops = TRUE, normalized = FALSE)
outDegreeU <- degree(graph_unclassified, mode="out", loops = TRUE, normalized = FALSE)
totalDegreeU <- degree(graph_unclassified)
inClosenessU <- closeness(graph_unclassified, mode='in')
outClosenessU <- closeness(graph_unclassified, mode='out')
totalClosenessU <- closeness(graph_unclassified)
betweennessU <- betweenness(graph_unclassified, directed = T)
eigenU <- evcent(graph_unclassified)
bonU <- bonpow(graph_unclassified)

sumU <- data.frame(inDegreeU, outDegreeU, totalDegreeU, inClosenessU, outClosenessU, totalClosenessU, betweennessU, eigenU)
sumU = sumU[,c(1:9)]

# correlate the measures
corrC <- cor(sumC)
corrU <- cor(sumU)