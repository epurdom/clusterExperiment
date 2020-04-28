data(simData)
res_km <- mainClustering(clusterFunction="kmeans", inputMatrix=simData,
               inputType="X", clusterArgs=list(k=3))
table(res_km, trueCluster)

res_mbkm <- mainClustering(clusterFunction="mbkmeans", inputMatrix=simData,
                         inputType="X", clusterArgs=list(k=3))
table(res_mbkm, trueCluster)

res_snn <- mainClustering(clusterFunction="snn", inputMatrix=simData,
                          inputType="X", clusterArgs=list(k=10))
table(res_snn, trueCluster)

### real data
set.seed(14456) ## for reproducibility, just in case
library(clusterExperiment)
data(fluidigmData) ## list of the two datasets (tophat_counts and rsem_tpm)
data(fluidigmColData)
se<-SummarizedExperiment(fluidigmData,colData=fluidigmColData)

system.time(ce<-clusterMany(se,clusterFunction=c("mbkmeans", "kmeans"),ks=5:10, minSizes=5,
                isCount=TRUE,reduceMethod="PCA",nReducedDims=c(5,15,50),run=TRUE))

system.time(cons1 <- makeConsensus(ce, proportion = 0.7))
system.time(cons2 <- makeConsensus(ce, proportion = 0.7, clusterFunction = "snn",
                                   clusterArgs = list(k=10)))

system.time(rsec1<-RSEC(se, isCount=TRUE, reduceMethod="PCA", nReducedDims=c(50,10), k0s=4:15,
              alphas=c(0.1,0.2,0.3),betas=c(0.8,0.9), minSizes=c(1,5), clusterFunction="hierarchical01",
              consensusProportion=0.7, consensusMinSize=5,
              dendroReduce="PCA",dendroNDims=50,
              mergeMethod="adjP",mergeCutoff=0.05,
))
