# Using KNN to solve NxN problem

## Background

Currently the results of subsampling (or consensus between clusters) is stored as a NxN matrix, with entries $s_{ij}$ being the number of times cell $i$ and $j$ are in the same cluster over the $B$ clusterings of which they both were assigned. 

*Note*: For subsampling, they may not be assigned because they were not both subsampled; for consensus they may not be assigned because one of them was a -1. 

Then `hier01` clusters this NxN matrix by taking entries $d_{ij}=1-s_{ij}$ and performing hierarchical clustering on the NxN distance matrix $D$. Afterward, we determine clusters from the hierarchical clustering by starting at the top of the tree and going down the tree until the resulting clustering satisfies the criteria that the cells in the cluster are sufficiently similar to each other, with a parameter $\alpha$ that controls the required similarity. There are two choices

1) $\forall i,j \in \mathcal{C}, s_{ij} \geq 1-\alpha$, i.e. the maximal allowable $d_{ij}$ for any $i,j$ in the cluster is $\alpha$. 
2) $\forall i \in \mathcal{C}, $mean(s_{ij}) \geq 1-\alpha$

Instead of calculating and storing the $NxN$ matrix, we would like to store the $NxB$ matrix of the cluster assignments. 

Davide proposed using KNN functionality to determine the neighbors of a cell $i$ whose distance to other cells is such that $d_{ij}<\alpha$. Then we would need to figure out how to get a cluster assignment out of that.  

## Proposal

1) We use `findNeighbors` from the `BiocNeighbors` package using `threshold=alpha`. 

**Problem** this function only works with Manhattan distance and Euclidean distance. We would need to see if the algorithm is even valid for our arbitrary distance (basically a Hamming distance) -- Aaron would probably know. Then we would need get Aaron to allow us to do a pull request and write C code for our distance function to be added to `BiocNeighbors`. Luckily it has already been updated to allow for two, so adding a third would hopefully not be a problem.

2) Convert these results into an object for `igraph`

3) Determine the clusters. There are three logical choices, that would not require any further parameters to be passed.
  a) Take all maximally connected components with the function `components` in `igraph`. This is a far looser requirement than any we have now.
  b) Take all maximal cliques (i.e. fully connected components) with the function `igraph_maximal_cliques`. This would basically correspond to the maximal option we currently have
  c) From each maximally connected component iteratively prune off vertices, based on the criteria that the remaining cluster must have all vertices with degree $< |V|/2$, where $|V|$ is the size of the cluster -- i.e. require them to be connected to at least half the other vertices in the cluster. This is similar to the mean criteria above, if we consider the median instead of the mean. (It would have to be iterative since removing a node will change the degree of all remaining nodes, though it could be made simpler by removing all such vertices at each iteration, though this is not exactly equivalent). However, note that in `hier01` our iterations go down the tree and split the existing cluster rather than pruning away cell by cell. 
  d) apply some community clustering algorithm to the $\alpha$-graph, which is equivalent to applying it to each maximally connected component for any algorithm worth anything. 
  
  
  