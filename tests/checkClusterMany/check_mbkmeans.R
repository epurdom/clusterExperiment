data(simData)
res_km <- mainClustering(clusterFunction="kmeans", inputMatrix=simData,
               inputType="X", clusterArgs=list(k=3))
table(res_km, trueCluster)

res_mbkm <- mainClustering(clusterFunction="mbkmeans",
                           inputMatrix=simData,
                           inputType="X", clusterArgs=list(k=3))
table(res_mbkm, trueCluster)

### real data
set.seed(14456) ## for reproducibility, just in case
#library(clusterExperiment)
data(fluidigmData) ## list of the two datasets (tophat_counts and rsem_tpm)
data(fluidigmColData)
se<-SummarizedExperiment(fluidigmData,colData=fluidigmColData)

system.time(ce<-clusterMany(se,clusterFunction=c("mbkmeans", "kmeans"),ks=5:10, minSizes=5,
                isCount=TRUE,reduceMethod="PCA",nReducedDims=c(5,15,50),run=TRUE))

system.time(cons1 <- makeConsensus(ce, proportion = 0.7))
system.time(cons2 <- makeConsensus(ce, proportion = 0.7,
                                   clusterFunction = "snn"))

system.time(rsec1<-RSEC(se, isCount=TRUE, reduceMethod="PCA", nReducedDims=c(50), k0s=4:7,
              alphas=c(0.1),betas=c(0.8), minSizes=c(5), clusterFunction="hierarchical01",
              consensusProportion=0.7, consensusMinSize=5,
              dendroReduce="PCA",dendroNDims=50,
              mergeMethod="none",
              subsampleArgs = list(clusterFunction="kmeans",
                                   resamp.n=100, samp.p=0.7)
))

system.time(rsec2<-RSEC(se, isCount=TRUE, reduceMethod="PCA", nReducedDims=c(50), k0s=4:7,
                        alphas=c(0.1),betas=c(0.8), minSizes=c(5), clusterFunction="snn", mainClusterArgs = list(algorithm = "louvain"),
                        consensusProportion=0.7, consensusMinSize=5,
                        dendroReduce="PCA",dendroNDims=50,
                        mergeMethod="none",
                        subsampleArgs = list(clusterFunction="kmeans",
                                             resamp.n=100, samp.p=0.7),
                        consensusArgs = list(clusterFunction="snn")
))
table(primaryCluster(rsec1), primaryCluster(rsec2))


tmp1 <- clusterSingle(simData, subsample = TRUE, sequential = FALSE,
              mainClusterArgs = list(clusterFunction = "hierarchical01",
                                     clusterArgs = list(alpha = 0.9)),
              subsampleArgs = list(clusterFunction = "kmeans",
                                   clusterArgs = list(k = 5),
                                   samp.p = 0.7,
                                   resamp.num = 100))

tmp2 <- clusterSingle(se[1:100,], subsample = TRUE, sequential = FALSE,
                      mainClusterArgs = list(clusterFunction = "snn",
                                             clusterArgs = list(alpha = 0.9)),
                      subsampleArgs = list(clusterFunction = "kmeans",
                                           clusterArgs = list(k = 5),
                                           samp.p = 0.7,
                                           resamp.num = 100))


tmp = clusterSingle(simData, subsample=TRUE, sequential=FALSE,
                    mainClusterArgs = list(
                      clusterFunction="hierarchical01",
                      clusterArgs=list(alpha=0.9)
                    ),
                    subsampleArgs = list(
                      clusterFunction = "kmeans",
                      clusterArgs = list(k = 5),
                      samp.p = 0.7,
                      resamp.num = 100))

tmp = clusterSingle(simData, subsample=TRUE, sequential=FALSE,
                    mainClusterArgs = list(
                      clusterFunction="snn",
                      clusterArgs=list(alpha=0.9)
                    ),
                    subsampleArgs = list(
                      clusterFunction = "kmeans",
                      clusterArgs = list(k = 5),
                      samp.p = 0.7,
                      resamp.num = 100))

tmp = clusterSingle(simData,
                    subsample=TRUE, sequential=FALSE,
                    mainClusterArgs = list(
                      clusterFunction="snn",
                      clusterArgs=list(alpha=0.2, algorithm="louvain")
                    ),
                    subsampleArgs = list(
                      clusterFunction = "kmeans",
                      clusterArgs = list(k = 5),
                      samp.p = 0.7,
                      resamp.num = 100))
