library(igraph)

#  Import Lazega friend adjacency matrix
lazega_friends_matrix <- as.matrix(read.csv(file.choose() ,header=TRUE,row.names=NULL,check.names=FALSE))

#  Create igraph object from this matrix, "undirected" will put a maximum of 1
#  for both edges (A->B and B->A) if either edge has a 1; other options are 
#  available, including "directed", "max", "min", etc 
lazega_friends_graph <- graph.adjacency(lazega_friends_matrix,mode="directed",weighted=NULL)

# import Lazega vertex attribute data (as a data.frame)
lazega_attributes <- read.csv(file.choose(), header=TRUE) # see  Lazega-atts.csv  on Courseworks

# attach the attributes to the Vertices of the igraph
vertex_attr(lazega_friends_graph, index=lazega_attributes$ID) <- lazega_attributes
# `index=lazega_attributes$ID" is only necessary if the order of the vertices and
#  the order of the attributes might be different. Better safe than sorry!

# calculate degree
lazega_attributes <- merge(lazega_attributes,  #  Merge lazega_attributes
                           data.frame(  # With a new data.frame
                             ID=V(lazega_friends_graph)$ID,  # Where the ID is the ID of each vertex
                             degree= degree(lazega_friends_graph) # and the degee is its degree
                           ),
                           by='ID')

head(lazega_attributes)

# calculate other centrality measures
lazega_attributes <- merge(lazega_attributes,  #  Merge lazega_attributes
                           data.frame(  # With a new data.frame
                             ID=V(lazega_friends_graph)$ID,  # Where the ID is the ID of each vertex
                             in.deg= degree(lazega_friends_graph, mode = c("in"), loops = TRUE, normalized = FALSE),
                             out.deg= degree(lazega_friends_graph, mode = c("out"), loops = TRUE, normalized = FALSE),
                             btwn= betweenness(lazega_friends_graph, directed = F),
                             close = closeness(lazega_friends_graph, mode = c("all")),
                             eigen <- evcent(lazega_friends_graph),
                             bon <- bonpow(lazega_friends_graph)
                           ),
                           by='ID')

names(lazega_attributes )

# clean up the dataframe
lazega_attributes = lazega_attributes[,c(1:14, 36)]

# change the bon centrality name
names(lazega_attributes)[names(lazega_attributes)=="bon....bonpow.lazega_friends_graph."] <- "bon"

# correlate the measures
cor(lazega_attributes)

# regress attributes on centrality measures
summary(lm(in.deg ~ status + out.deg + gender + practice, lazega_attributes))

## make a pretty graph of your network

## start the graph ##
set.seed(12)
l <- layout.kamada.kawai(lazega_friends_graph)

# Plot undecorated first.
par(mfrow=c(1,1))
oldMargins<-par("mar")
par(mar=c(1,1,1,1))
### par(mar=oldMargins) ### to return to default ...

# Differentiate two gender by shape
V(lazega_friends_graph)[V(lazega_friends_graph)$gender== 1]$shape <- "circle"
V(lazega_friends_graph)[V(lazega_friends_graph)$gender== 2]$shape <- "rectangle"

# Differentiate two status by color.
V(lazega_friends_graph)[V(lazega_friends_graph)$status== 1]$color <- "white"
V(lazega_friends_graph)[V(lazega_friends_graph)$status== 2]$color <- "dodgerblue"

# Differentiate two practices by label color.
V(lazega_friends_graph)[V(lazega_friends_graph)$practice== 1]$label.color <- "red"
V(lazega_friends_graph)[V(lazega_friends_graph)$practice== 2]$label.color <- "black"

# Size node by in-degree.
V(lazega_friends_graph)$size <- 4*sqrt(degree(lazega_friends_graph, mode="in"))
V(lazega_friends_graph)$size2 <- V(lazega_friends_graph)$size * .5

# Size of node label by out-degree.
V(lazega_friends_graph)$label.cex <- 2.5 * degree(lazega_friends_graph, mode="out") / max(degree(lazega_friends_graph, mode="out"))

# Shrink arrows
plot(lazega_friends_graph, layout=l, edge.arrow.size=.3)

## to suppress the edges, add edge.color="white" ##
## plot(testmg, layout=l, edge.color="white")