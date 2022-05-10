##Code for timing the results of the tests:
# testOutput<-testOutput<-test_dir("testthat/",reporter=ListReporter)
# timingDf<-do.call("rbind", lapply(testOutput, function(x){data.frame(file=x$file, test=x$test, time=round(x$real,2))}))
# timingDf<-timingDf[order(timingDf$time),]
# write.table(file="testTimings.txt", x=timingDf)

###Note: any changes to this file should be at the END so as to not mess up the seed calls.
library(clusterExperiment)
# library(devtools)
# load_all()
data(simData, envir = environment())
if(ncol(simData) != 300) {
  stop("not current version of simData")
  #get all kinds of annoyances because using old version.
  #Can delete this once package is stabilized.
}

#################################
###Simple, trivial sized objects for testing:
# mat=20x15 matrix of data with row and column names
# labMat 15 x 2 matrix of numberic clusters, labeled 2-5, plus -1,-2 values with column names
# sData 15 x 3 matrix of data on the samples, a factor, continuous, and character variable with column names
# gData 20 x 4 matrix of data on the features: factor, continuous, character (like sData) with column names
# mData list of length 2 of random information to be metaData
# se a Summarized experiment with rowData, colData, and metaData slots
#################################
set.seed(54)
mat <- matrix(data=rnorm(20*15), ncol=15)
mat[1,1]<- -1 #force a negative value
colnames(mat)<-paste("Sample",1:ncol(mat))
rownames(mat)<-paste("Gene",1:nrow(mat))
set.seed(2325)
se <- simSEObject(mat)
gData<-rowData(se)
#dissMat<-as.matrix(dist(t(mat)))

# Create matrix of clusters:
numLabels <- as.character(gl(5, 3))
numLabels[c(1:2)]<- c("-1","-2") #make sure some not assigned
numLabels<-factor(numLabels)
chLabels<-rep(LETTERS[1:5],each=3)
chLabels[c(2:3)]<- c("-1","-2") #make sure some not assigned
labMat<-cbind(as.numeric(as.character(numLabels)),as.numeric(as.character(numLabels)))
colnames(labMat)<-c("Cluster1","Cluster2")

# Construct CE Objects
cc <- ClusterExperiment(mat, labMat, transformation = function(x){x})
ccSE<-ClusterExperiment(se, labMat, transformation=function(x){x})


#################################
###Larger sized objects based on simData/simCount:
#################################
data(simData)
if(ncol(simData)!=300) stop("not current version of simData") #got all kinds of annoyances from using old version.
set.seed(2415)
seSimData <-simSEObject(simData)
# Same meta data, but with the counts
set.seed(2415)
seSimCount <-simSEObject(simData)
data(simCEData)

#################################
###small object based on simData/simCount (same size as trivial data)
### 15 samples from each of groups (including -2,-1)
###
#################################
set.seed(53289)
whSamp<-unlist(tapply(1:nSamples(ceSimCount),primaryCluster(ceSimCount),function(x){sample(x,size=3)})) #15
smSimData<-simData[1:20,whSamp]
smSimCount<-simCount[1:20,whSamp]
smSimCE<-ceSimCount[1:20,whSamp]
smSimSE <- seSimData[1:20,whSamp]



#################################
###Make reduce dimensions and filters
###... needed for creation of hdf5 object AND testing reducedDims
#################################
sce<-as(se,"SingleCellExperiment")
sceFull<-sce
clusterExperiment:::filterStats(sceFull,type=c("Filter1","Filter2"))<-matrix(rnorm(2*nrow(sce)),ncol=2)
reducedDim(sceFull,type="Red1")<-matrix(rnorm(2*ncol(sce)),ncol=2)


sceSimData<-as(seSimData,"SingleCellExperiment")
sceSimDataDimRed<-sceSimData
pca_data <- prcomp(t(assay(sceSimData)),scale=TRUE,center=TRUE)
tsne_data <- matrix(rnorm(NCOL(sceSimData)*2),ncol=2)
reducedDims(sceSimDataDimRed) <- SimpleList(PCA=pca_data$x, TSNE=tsne_data)
clusterExperiment:::filterStats(sceSimDataDimRed,type=c("Filter1","Filter2"))<-matrix(rnorm(2*nrow(sceSimDataDimRed)),ncol=2)

#####################
## Create hdf5 SCE version
## Note is matrix of doubles....
#####################
# hdfSCE<-HDF5Array::saveHDF5SummarizedExperiment(sceSimDataDimRed, dir="sceRedDem.h5", replace=TRUE)
# hdfObj<-HDF5Array::saveHDF5SummarizedExperiment(sceSimData, dir="sce.h5", replace=TRUE)

