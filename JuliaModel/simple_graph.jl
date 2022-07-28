using SimpleWeightedGraphs, Graphs

# First we construct a weighted, directed graph
g_weighted = SimpleWeightedDiGraph(G)

# For later use we extract the edge.weight attributes
# . is the broadcasting operator and gets the attribute :weight for every edge
edge_weights = getfield.(collect(edges(g_weighted)), :weight)

# we promote the g_weighted graph as a directed graph (weights of the edges are included in parameters)
g_directed = SimpleDiGraph(g_weighted)