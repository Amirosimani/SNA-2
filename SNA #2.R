library(igraph)

### importing data ----
#  Import adjacency matrix
hillary_matrix <- as.matrix(read.csv(file.choose() ,header=TRUE,row.names=NULL,check.names=FALSE))
hillary_matrix <- hillary_matrix[,-1]

#  Create igraph object from this matrix
hillary_graph <- graph.adjacency(hillary_matrix,mode="directed", weighted = TRUE)

###### attributes

### overall topography of the network ----
vcount(hillary_graph)
ecount(hillary_graph)

g_density = graph.density(hillary_graph, loops=T)


### centrality measures ----
inDegree <- degree(hillary_graph, mode="in", loops = TRUE, normalized = FALSE)
outDegree <- degree(hillary_graph, mode="out", loops = TRUE, normalized = FALSE)
totalDegree <- degree(hillary_graph)
inCloseness <- closeness(hillary_graph, mode='in')
outCloseness <- closeness(hillary_graph, mode='out')
totalCloseness <- closeness(hillary_graph)
betweenness <- betweenness(hillary_graph, directed = T)
eigen <- evcent(hillary_graph)
bon <- bonpow(hillary_graph)

forum <- data.frame(inDegree, outDegree, totalDegree, inCloseness, outCloseness, totalCloseness, betweenness, bon, eigen)

# clean up the dataframe
forum = forum[,c(1:9)]

# correlate the measures
cor <- cor(forum)
