# Using KNN to solve NxN problem

## Background

Currently the results of subsampling (or consensus between clusters) is stored as a NxN matrix, with entries $s_{ij}$ being the proportion of times cell $i$ and $j$ are in the same cluster over the $B$ clusterings of which they both were assigned. 

*Note*: For subsampling, they may not be assigned because they were not both subsampled; for consensus they may not be assigned because one of them was a -1. 

Then `hier01` clusters this NxN matrix by taking entries $d_{ij}=1-s_{ij}$ and performing hierarchical clustering on the NxN distance matrix $D$. Afterward, we determine clusters from the hierarchical clustering by starting at the top of the tree and going down the tree until the resulting clustering satisfies the criteria that the cells in the cluster are sufficiently similar to each other, with a parameter $\alpha$ that controls the required similarity. There are two choices

1) $\forall i,j \in \mathcal{C}, s_{ij} \geq 1-\alpha$, i.e. the maximal allowable $d_{ij}$ for any $i,j$ in the cluster is $\alpha$. 

2) $\forall i \in \mathcal{C}, mean(s_{ij}) \geq 1-\alpha$ [default]

Instead of calculating and storing the $NxN$ matrix, we would like to store the $NxB$ matrix of the cluster assignments. 

Davide proposed using KNN functionality to determine the neighbors of a cell $i$ whose distance to other cells is such that $d_{ij}<\alpha$. Then we would need to figure out how to get a cluster assignment out of that.  

## Proposal

1) We use `findNeighbors` from the `BiocNeighbors` package using `threshold=alpha`. 

**To Do** This function only works with Manhattan distance and Euclidean distance. We would need to see if the algorithm is even valid for our arbitrary distance (basically a Hamming distance) -- Aaron would probably know. Then we would need get Aaron to allow us to do a pull request and write C code for our distance function to be added to `BiocNeighbors`. Luckily it has already been updated to allow for two, so adding a third would hopefully not be a problem.

2) Convert these results into an object for `igraph`. 

**To Do** Not clear that there is an existing function to do this. Might be relevant add on to `BiocNeighbors`. More generally if we create a wrapper package for the clustering methods that we were talking about, I think interface between igraph and `BiocNeighbors` may be needed, but I haven't looked into this. 

3) Determine the clusters. There are several logical choices, only some of which would not require any further parameters to be passed.

 a) Take all maximally connected components with the function `components` in `igraph`. This is a far looser requirement than any we have now.
  
  b) Take all maximal cliques (i.e. fully connected components) with the function `igraph_maximal_cliques`. This would basically correspond to the maximal option we currently have
  
  c) From each maximally connected component iteratively prune off vertices, based on the criteria that the remaining cluster must have all vertices with degree $< |V|/2$, where $|V|$ is the size of the cluster -- i.e. require them to be connected to at least half the other vertices in the cluster. This is similar to the mean criteria above, if we had considered the median instead of the mean. (It would have to be iterative since removing a node will change the degree of all remaining nodes). Could use some varient of this without iteration, where we look at degree distribution per connected components and pick some value for degree edge; perhaps by looking at distribution it's possible to find a degree number that would guaranteed mathematically that the resulting vertices have degree at least $<|V|/2$ without iteration (e.g. if removing 25% of vertices means that the remaining had degree $>.75|V|$ on full graph, then clearly pruned graph would have degree $>|V|/2$). 
  
  However, note that in `hier01` our iterations go down the tree and split the existing cluster rather than pruning away cell by cell. Without that difference, this pruning could easily devolve into b) with much worse speed and perhaps even worse (after all if a connected component consists of 2 maximal cliques, then maximal cliques would find both, while iterative would only get one of them). So seems like its important to have ability to "capture" those that are discarded, and not clear how to do that.
  
  d) apply some community clustering algorithm to the $\alpha$-graph, which is equivalent to applying it to each maximally connected component (for any algorithm worth anything). This would likely require some parameter to the clustering. Maybe question is how bad this parameter is relative to our current (arbitrary) options; advantage of our current method is at least based on a guarantee on resulting clustering regarding $\alpha$ in some way, which a clustering wouldn't necessarily. 
  
It is unclear how different a) and b) will from each other in practice -- i.e. how much extra is added by including in any cell with small alpha to at least one cell. Clearly if not, none of these differences matter. But seems like we need some kind of method like our current one that is somewhere in between the two, since we don't know.
  
## Updating findNeighbors

### General outline of the function
`findNeighbors` is a upper-level function that finds all neighbors within a given distance of each point. It takes a parameter `BiocNeighborParam` to define the algorithm to do so (this is a argument of class `BiocNeighborParam`). `findKNN` is the corresponding function that finds the k-nearest neighbors of each point. 

`findNeighbors` has two possible options for the algorithm, which is basically a choice of how to most efficiently organize the data in the search: "Kmknn" or "VPTree". `findNeighbors` essentially dispatches `rangeFindKmknn` or `rangeFindVptree`. Note that these functions allow the `threshold` argument to be a vector, specifying a different threshold for each point, though not clear this is useful for us. `findNeighbors` will return not just the indices of the points, but also their distances (if they are less than the threshold).  (`findKNN` also has approximate methods `Annoy` and `Hnsw` but these do not appear to be options for finding within a particular distance; it would be important to consider in comparison methods, like Seurat that build a knn neighbor graph, whether the speed increases they have for large datasets are based on using approximations or not). 

(I would note however, that the  help function for `BiocNeighborParam` only lists "2 concrete subclasses",  `KmknnParam` and `AnnoyParam`. It's not clear why the help doesn't list the other two options as subclasses.)

biocNeighbors has two methods for preprocessing data for fast (exact) computation (from :

* KMKNN (k-means nearest neighbor): In the KMKNN algorithm (Wang, 2012), k-means clustering is first applied to the data points usingthe square root of the number of points as the number of cluster centers. The cluster assignment anddistance to the assigned cluster center for each point represent the KMKNN indexing information.This speeds up the nearest neighbor search by exploiting the triangle inequality between clustercenters, the query point and each point in the cluster to narrow the search space. The advantage ofthe KMKNN approach is its simplicity and minimal overhead, resulting in performance improve-ments over conventional tree-based methods for high-dimensional data where most points need tobe searched anyway.  It is also trivially extended to find all neighbors within a threshold distancefrom a query point.  Note that KMKNN operates much more naturally with Euclidean distances, so your mileage may vary when using it with Manhattan distances.
	- The pre-step runs kmeans on the data (`buildKmknn` using `stat::kmeans`, all R code)
	- Then `.range_find_kmknn` is run, which is a call to C code: `.Call(cxx_range_find_kmknn,...)`
* Vantage-point tree (VPTree): In a VP tree (Yianilos, 1993), each node contains a subset of points and has a defined threshold distance (usually the median distance to all points in the subset). The left child of this node containsthe further subset of points within the radius, while the right child contains the remaining points inthe subset that are outside.  The nearest neighbor search traverses the tree and exploits the triangle inequality between query points, node centers and thresholds to narrow the search space. VP trees are often faster than more conventional KD-trees or ball trees as the former uses the points them-selves as the nodes of the tree, avoiding the need to create many intermediate nodes and reducing the total number of distance calculations.  Like KMKNN, it is also trivially extended to find all neighbors within a threshold distance from a query point. "A vantage-point tree (or VP tree) is a metric tree that segregates data in a metric space by choosing a position in the space (the "vantage point") and partitioning the data points into two parts: those points that are nearer to the vantage point than a threshold, and those points that are not. By recursively applying this procedure to partition the data into smaller and smaller sets, a tree data structure is created where neighbors in the tree are likely to be neighbors in the space... The vantage-point tree is particularly useful in dividing data in a non-standard metric space into a metric tree... All it needs is the distance function that satisfies the properties of the metric space... The time cost to build a Vantage-Point tree is approximately O(n log n)" (Wikipedia). 
	- The pre-step builds the VP Tree (`buildVptree` which is a C call, `.Call(cxx_build_vptree,...)`
	- Then call to `.range_find_vptree` (`.Call(cxx_range_find_vptree,...)`)

### Our distance

We need to be sure that our proposed distance satifies the triangle-inequality / defines a proper metric. Otherwise, neither of these algorithms may work... Hamming distance is a proper metric, but we are doing an alternative of it, to deal with NA values (i.e. not use the NA values, but then normalize the distance to account for the number of positions shared).

*How important is this?* For subsampling, we can have samples that were not subsampled togeter, and thus have NA for their clustering. Our default implementation doesn't actually do this -- instead it assignes every sample, not only the subsampled, based on the clusters defined by the subsampled samples. But if user picked the other option, then it would be a problem. 

For consensus clustering, we get -1 values, and we do not count that against the sample. Namely, we ignore the -1, get the proportion of times together; we then filter out those with too many -1. 

If dealing with the proportion becomes a problem, rather than straightforward Hamming distance, then we could maybe reconsider these options. 

*Algorithms* Clearly the kmknn is not a good option for our distance, since kmeans on our $N x B$ matrix is not going to be terribly informative. 