library(clusterExperiment)
data(simData)
if(ncol(simData) != 300) {
  stop("not current version of simData")
  #get all kinds of annoyances because using old version.
  #Can delete this once package is stabilized.
}
options(getClass.msg=FALSE) #get rid of annoying messages about cache so not printed on build

## make sure the tests are reproducible
set.seed(23)

#################################
###Simple, trivial sized objects for testing:
# mat=20x15 matrix of data with row and column names
# labMat 15 x 2 matrix of numberic clusters, labeled 2-5, plus -1,-2 values with column names
# sData 15 x 3 matrix of data on the samples, a factor, continuous, and character variable with column names
# gData 20 x 4 matrix of data on the features: factor, continuous, character (like sData) with column names
# mData list of length 2 of random information to be metaData
# se a Summarized experiment with rowData, colData, and metaData slots
#################################
mat <- matrix(data=rnorm(20*15), ncol=15)
mat[1,1]<- -1 #force a negative value
colnames(mat)<-paste("Sample",1:ncol(mat))
rownames(mat)<-paste("Gene",1:nrow(mat))
numLabels <- as.character(gl(5, 3))
numLabels[c(1:2)]<- c("-1","-2") #make sure some not assigned
numLabels<-factor(numLabels)
chLabels<-rep(LETTERS[1:5],each=3)
chLabels[c(2:3)]<- c("-1","-2") #make sure some not assigned
labMat<-cbind(as.numeric(as.character(numLabels)),as.numeric(as.character(numLabels)))
colnames(labMat)<-c("Cluster1","Cluster2")
sData<-data.frame(sample(letters[2:5],size=NCOL(mat),replace=TRUE),sample(2:5,size=NCOL(mat),replace=TRUE))
sData<-data.frame(sData,sample(LETTERS[2:5],size=NCOL(mat),replace=TRUE),stringsAsFactors=FALSE)
gData<-data.frame(sample(letters[2:5],size=NROW(mat),replace=TRUE),sample(2:5,size=NROW(mat),replace=TRUE))
gData<-data.frame(gData,sample(LETTERS[2:5],size=NROW(mat),replace=TRUE),stringsAsFactors=FALSE)
colnames(sData)<-c("A","B","C")
colnames(gData)<-c("a","b","c")
mData<-list(first=c(1,2,3),second=c("Information"))
se <- SummarizedExperiment(mat,colData=sData,rowData=gData,metadata=mData)
cc <- clusterExperiment(mat, labMat, transformation = function(x){x})
ccSE<-clusterExperiment(se,labMat,transformation=function(x){x})

#################################
###Larger sized objects based on simData/simCount:
#################################
data(simData)
if(ncol(simData)!=300) stop("not current version of simData") #get all kinds of annoyances because using old version.
simSData<-data.frame(sample(letters[2:5],size=NCOL(simData),replace=TRUE),sample(2:5,size=NCOL(simData),replace=TRUE))
simSData<-data.frame(simSData,sample(LETTERS[2:5],size=NCOL(simData),replace=TRUE),stringsAsFactors=FALSE)
colnames(simSData)<-c("A","B","C")
gSimData<-data.frame(sample(letters[2:5],size=NROW(simData),replace=TRUE),sample(2:5,size=NROW(simData),replace=TRUE))
gSimData<-data.frame(gSimData,sample(LETTERS[2:5],size=NROW(simData),replace=TRUE),stringsAsFactors=FALSE)
colnames(gSimData)<-c("a","b","c")

seSimData <- SummarizedExperiment(simData,colData=simSData,rowData=gSimData,metadata=mData)
seSimCount <- SummarizedExperiment(simCount,colData=simSData,rowData=gSimData,metadata=mData)

test<- clusterMany(simCount,dimReduce="PCA",nPCADims=c(5,10,50), isCount=TRUE,
                         clusterFunction="pam",ks=2:4,findBestK=c(TRUE,FALSE))
test<-addClusters(test,sample(2:5,size=NCOL(simData),replace=TRUE),clusterTypes="User")
clMatNew<-apply(clusterMatrix(test),2,function(x){
    wh<-sample(1:nSamples(test),size=10)
    x[wh]<- -1
    wh<-sample(1:nSamples(test),size=10)
    x[wh]<- -2
    return(x)
})
#make a new object with -1 values
ceSim<-clusterExperiment(seSimCount,clMatNew,transformation=function(x){log2(x+1)})
clusterTypes(ceSim)<-clusterTypes(test)
rm(test)
#################################
###small object based on simData/simCount (same size as trivial data)
### 15 samples from each of groups (including -2,-1)
###
#################################
whSamp<-unlist(tapply(1:nSamples(ceSim),primaryCluster(ceSim),function(x){sample(x,size=3)})) #15
smSimData<-simData[1:20,whSamp]
smSimCount<-simCount[1:20,whSamp]
smSimCE<-ceSim[1:20,whSamp]
smSimSE <- seSimData[1:20,whSamp]
