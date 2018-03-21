library(scran)
library(zinbwave)
library(Rtsne)
library(igraph)
library(profvis)
library(devtools)
library(testthat)
library(pryr)
library(microbenchmark)

load("tests/checkClusterMany/combined_filtered_all_batches.rda")
load("tests/checkClusterMany/combined_zinbwave_all_batches.rda")

W <- getW(zinb)
rownames(W) <- colnames(filtered)
colnames(W) <- paste0("W", 1:10)
WW <- t(W)

## profile
source("tests/checkClusterMany/subsample_for_profile.R")
profvis(subsample(getBuiltInFunction("kmeans"), x = WW,
                  clusterArgs = list(k=5),
                  samp.p = 0.7, resamp.num = 10,
                  largeDataset = TRUE, whichImplementation = "Cmemory"))

# trap_wrap <- function(x, k, steps = 4, ...) {
#   snn <- buildSNNGraph(x, k = k, d = NA, transposed = FALSE)
#   res <- cluster_walktrap(snn, steps = steps)
#   return(res$membership)
# }
#
# internalFunctionCheck(trap_wrap, inputType = "X", algorithmType = "K",
#                       outputType="vector")
#
# SNN <- ClusterFunction(trap_wrap, inputType = "X", algorithmType = "K",
#                        outputType="vector")

library(clusterExperiment)

microbenchmark(
  original={clusterSingle(t(W), subsample = TRUE, sequential = FALSE,
                          mainClusterArgs = list(clusterFunction = "hierarchical01",
                                                 clusterArgs = list(alpha = 0.9)),
                          subsampleArgs = list(clusterFunction = "kmeans",
                                               clusterArgs = list(k = 5),
                                               samp.p = 0.7,
                                               resamp.num = 100))},
  large_r={clusterSingle(t(W), subsample = TRUE, sequential = FALSE,
                         mainClusterArgs = list(clusterFunction = "hierarchical01",
                                                clusterArgs = list(alpha = 0.9)),
                         subsampleArgs = list(clusterFunction = "kmeans",
                                              clusterArgs = list(k = 5),
                                              samp.p = 0.7,
                                              resamp.num = 100,
                                              largeDataset = TRUE,
                                              whichImplementation = "R"))},
  large_c1={clusterSingle(t(W), subsample = TRUE, sequential = FALSE,
                         mainClusterArgs = list(clusterFunction = "hierarchical01",
                                                clusterArgs = list(alpha = 0.9)),
                         subsampleArgs = list(clusterFunction = "kmeans",
                                              clusterArgs = list(k = 5),
                                              samp.p = 0.7,
                                              resamp.num = 100,
                                              largeDataset = TRUE,
                                              whichImplementation = "Csimple"))},
  large_c2={clusterSingle(t(W), subsample = TRUE, sequential = FALSE,
                         mainClusterArgs = list(clusterFunction = "hierarchical01",
                                                clusterArgs = list(alpha = 0.9)),
                         subsampleArgs = list(clusterFunction = "kmeans",
                                              clusterArgs = list(k = 5),
                                              samp.p = 0.7,
                                              resamp.num = 100,
                                              largeDataset = TRUE,
                                              whichImplementation = "Cmemory"))},
  times = 1L
)








## first let's make sure that we get the same results, then we can memory proof them

devtools::load_all()
set.seed(123)
system.time(lds1 <- clusterSingle(t(W), subsample = TRUE, sequential = FALSE,
                                 mainClusterArgs = list(clusterFunction = "hierarchical01",
                                                        clusterArgs = list(alpha = 0.9)),
                                 subsampleArgs = list(clusterFunction = "kmeans",
                                                      clusterArgs = list(k = 5),
                                                      samp.p = 0.1,
                                                      resamp.num = 10,
                                                      largeDataset = TRUE,
                                                      whichImplementation = "Csimple")))

set.seed(123)
system.time(lds2 <- clusterSingle(t(W), subsample = TRUE, sequential = FALSE,
                                    mainClusterArgs = list(clusterFunction = "hierarchical01",
                                                           clusterArgs = list(alpha = 0.9)),
                                    subsampleArgs = list(clusterFunction = "kmeans",
                                                         clusterArgs = list(k = 5),
                                                         samp.p = 0.1,
                                                         resamp.num = 10,
                                                         largeDataset = TRUE,
                                                         whichImplementation = "Cmemory")))

expect_equivalent(master, lds)

set.seed(123)
profvis(master <- clusterSingle(t(W), subsample = TRUE, sequential = FALSE,
                                    mainClusterArgs = list(clusterFunction = "hierarchical01",
                                                           clusterArgs = list(alpha = 0.9)),
                                    subsampleArgs = list(clusterFunction = "kmeans",
                                                         clusterArgs = list(k = 5),
                                                         samp.p = 0.7,
                                                         resamp.num = 100)))

set.seed(123)
profvis(lds <- clusterSingle(t(W), subsample = TRUE, sequential = FALSE,
                                 mainClusterArgs = list(clusterFunction = "hierarchical01",
                                                        clusterArgs = list(alpha = 0.9)),
                                 subsampleArgs = list(clusterFunction = "kmeans",
                                                      clusterArgs = list(k = 5),
                                                      samp.p = 0.7,
                                                      resamp.num = 100,
                                                      largeDataset = TRUE)))



classY <- sample(1:5, 20000, replace=TRUE)
profvis({
  D2 <- lapply(seq_len(length(classY)), function(i) classY == classY[i])
  D1 <- outer(classY, classY, function(a, b) a==b)
  D3 <- unlist(D2)
  D4 <- as.vector(D1)
})


master <- clusterSingle(t(W), subsample = TRUE, sequential = FALSE,
                        mainClusterArgs = list(clusterFunction = "hierarchical01",
                                               clusterArgs = list(alpha = 0.9)),
                        subsampleArgs = list(clusterFunction = "kmeans",
                                             clusterArgs = list(k = 5),
                                             samp.p = 0.1,
                                             resamp.num = 10, largeDataset=TRUE))


#############################


otherIds<-function(idx,clustVec,clustLeng){
  m<-which(clustVec==idx)
  if(length(m)>1) stop("ids clustered in more than one cluster")
  if(length(m)==0) return(NA) #sample not ever clustered
  if(length(m)==1){
    ends<-cumsum(clustLeng)
    begins<-cumsum(c(1,head(clustLeng,-1)))
    whCluster<-which(m<=ends & m>=begins)
    if(length(whCluster)>1 | length(whCluster)==0) stop("error in coding: finding range of clusterids")
    return(clustVec[seq(begins[whCluster],ends[whCluster],by=1)])
  }
}

otherIds2 <- function(idx, classX) {
  return(which(classX == classX[idx]))
}

otherIds3 <- Vectorize(otherIds2, vectorize.args = "idx", SIMPLIFY = FALSE)

searchForPairs<-function(ii,clusterList){

  #get list of those indices sample ii was sampled
  whHave<-which(sapply(clusterList,function(ll){ii%in%ll$clusterIds}))

  #calculate number of times sampled with (denominator)
  sampledWithTab<-table(unlist(lapply(clusterList[whHave],.subset2,"clusterIds")))

  #get those indices clustered with and tabulate
  clusterWith<-lapply(clusterList[whHave],function(ll){
    otherIds(idx=ii,clustVec=ll$clusterIds,clustLeng=ll$clusterLengths)
  })

  clusterWithTab<-table(unlist(clusterWith))
  jointNames<-as.character(1:N)
  whLower<-which(as.integer(as.numeric(jointNames))<ii)

  return(as.integer(clusterWithTab[jointNames][whLower])/as.integer(sampledWithTab[jointNames][whLower]))
  #old code that made Nx3 matrix and subset it:
  #             #make a N x 3 matrix summarizing the results for all indices
  # out<-cbind(idx=as.integer(as.numeric(jointNames)), together=as.integer(clusterWithTab[jointNames]), total=as.integer(sampledWithTab[jointNames]))
  # #keep only those in the lower triangle
  #             out<-out[out[,"idx"]<ii,,drop=FALSE]
  #             return(out[,"together"]/out[,"total"])

  #thoughts about alternative code if not need save NxN matrix...
  #jointNames<-names(sampledWithTab) #if manage to not save NxN matrix, could use this to return only those that actually present
  #out<-out[!is.na(out[,"together"]),,drop=FALSE] #if manage to not save NxN matrix, could use this to return only those that actually present; but then need to not return proportions, but something else.
}

searchForPairs2<-function(ii, clusterList, ncores = 1){

  #get list of those indices sample ii was sampled
  whHave<-which(sapply(clusterList,function(ll){ii%in%ll$clusterIds}))

  #calculate number of times sampled with (denominator)
  sampledWithTab<-table(unlist(lapply(clusterList[whHave],.subset2,"clusterIds")))

  #get those indices clustered with and tabulate
  if(ncores == 1) {
    clusterWith<-lapply(clusterList[whHave],function(ll){
      otherIds2(idx=ii, classX = ll$classX)
    })
  } else {
    clusterWith<-parallel::mclapply(clusterList[whHave],function(ll){
      otherIds2(idx=ii, classX = ll$classX, mc.cores = ncores)
    })
  }

  clusterWithTab<-table(unlist(clusterWith))
  jointNames<-as.character(1:N)
  whLower<-which(as.integer(as.numeric(jointNames))<ii)

  return(as.integer(clusterWithTab[jointNames][whLower])/as.integer(sampledWithTab[jointNames][whLower]))
  #old code that made Nx3 matrix and subset it:
  #             #make a N x 3 matrix summarizing the results for all indices
  # out<-cbind(idx=as.integer(as.numeric(jointNames)), together=as.integer(clusterWithTab[jointNames]), total=as.integer(sampledWithTab[jointNames]))
  # #keep only those in the lower triangle
  #             out<-out[out[,"idx"]<ii,,drop=FALSE]
  #             return(out[,"together"]/out[,"total"])

  #thoughts about alternative code if not need save NxN matrix...
  #jointNames<-names(sampledWithTab) #if manage to not save NxN matrix, could use this to return only those that actually present
  #out<-out[!is.na(out[,"together"]),,drop=FALSE] #if manage to not save NxN matrix, could use this to return only those that actually present; but then need to not return proportions, but something else.
}

searchForPairs3 <- Vectorize(searchForPairs2, vectorize.args = "ii", SIMPLIFY = FALSE)

## For each subsample (list) find the pairs
searchForPairs4 <- function(clusterMat) {

  retval <- matrix(ncol=N, nrow=N)

  for(i in seq(2, ncol(clusterMat))) {
    for(j in seq_len(i)) {
      s <- sum(clusterMat[,i] == clusterMat[,j], na.rm = TRUE)
      tot_na <- sum(!(is.na(clusterMat[,i]) | is.na(clusterMat[,j])))
      retval[i, j] <- s / tot_na
    }
  }

  return(retval)

}

library(testthat)
set.seed(112)
N <- 100
classX <- sample(1:10, N, replace = TRUE)
classY <- sample(1:10, N, replace = TRUE)
clusterIdList1 <- tapply(1:N, classX, identity, simplify=FALSE)
clusterIds1 <- unname(unlist(clusterIdList1))
clusterLengths1 <- as.integer(tapply(1:N, classX, length))
clusterIdList2 <- tapply(1:N, classY, identity, simplify=FALSE)
clusterIds2 <- unname(unlist(clusterIdList2))
clusterLengths2 <- as.integer(tapply(1:N, classY, length))

clusterList <- list(list(clusterIds=clusterIds1, clusterLengths=clusterLengths1, classX = classX),
              list(clusterIds=clusterIds2, clusterLengths=clusterLengths2, classX = classY))

expect_equal(lapply(1:N, otherIds2, classX), lapply(1:N, otherIds, clusterIds1, clusterLengths1))
expect_equal(lapply(1:N, otherIds, clusterIds1, clusterLengths1), otherIds3(1:N, classX))

expect_equal(lapply(2:N,function(jj){searchForPairs(jj,clusterList=clusterList)}),
             lapply(2:N,function(jj){searchForPairs2(jj,clusterList=clusterList)}))

expect_equal(lapply(2:N,function(jj){searchForPairs(jj,clusterList=clusterList)}), searchForPairs3(2:N, clusterList))

pairList <- searchForPairs3(2:N, clusterList)

#Create NxN matrix
Dbar<-matrix(0,N,N)
Dbar[upper.tri(Dbar, diag = FALSE)]<-unlist(pairList)
Dbar<-Dbar+t(Dbar)
Dbar[is.na(Dbar)]<-0
diag(Dbar)<-1
Dbar_lower <- Dbar[lower.tri(Dbar, diag=FALSE)]

clusterMat <- sapply(clusterList, function(x) x$classX)
Dbar2 <- search_pairs(t(clusterMat))
Dbar2_lower <- Dbar2[lower.tri(Dbar2, diag=FALSE)]

Dbar3 <- searchForPairs4(t(clusterMat))
Dbar3_lower <- Dbar3[lower.tri(Dbar3, diag=FALSE)]

expect_equal(Dbar_lower, Dbar3_lower)

library(microbenchmark)

microbenchmark(
  original={lapply(2:N,function(jj){searchForPairs(jj,clusterList=clusterList)})},
  vectorized={searchForPairs3(2:N, clusterList)},
  loopR={searchForPairs4(t(clusterMat))},
  loopCpp={search_pairs(t(clusterMat))}
)

# tapply(1:N, mat[,1], identity, simplify = FALSE)
library(cluster)

set.seed(1222)
N <- 25
x <- rbind(cbind(rnorm(10,0,0.5), rnorm(10,0,0.5)),
           cbind(rnorm(15,5,0.5), rnorm(15,5,0.5)))
mat <- t(sapply(2:4, function(k) pam(x, k=k)$clustering))
mat
mat2 <- mat[,1:3]
mat2[2,2] <- NA
mat2
mat3 <- mat2
mat3[2,3] <- NA
mat3
mat4 <- mat2
mat4[1,3] <- NA
mat4

search_pairs(mat2)
search_pairs(mat3)
search_pairs(mat4)

## TODO:
## 1. merge elizabeth branch
## 2. three options: old (R), my C++, Elizabeth's C++
##




