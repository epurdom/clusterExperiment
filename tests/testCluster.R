library(clusterCells)
data(simData)

#basic check over different dimensions and parameters -- without resampling
ps<-c(5,10,50)
names(ps)<-paste("npc=",ps,sep="")
pcaData<-stats::prcomp(simData, center=TRUE, scale=TRUE)
#check how many and what runs user choices will imply:
outPCADims <- compareChoices(lapply(ps,function(p){pcaData$x[,1:p]}), 
                              clusterMethod="pam",
                              ks=2:4,findBestK=c(TRUE,FALSE),removeSil=c(TRUE,FALSE),
                              run=FALSE)
#make names shorter for plotting
colnames(cl$clMat)<-gsub("TRUE","T",colnames(cl$clMat))
colnames(cl$clMat)<-gsub("FALSE","F",colnames(cl$clMat))
colnames(cl$clMat)<-gsub("k=NA,","",colnames(cl$clMat))
par(mar=c(2,10,1,1))
plotTracking(cl$clMat,axisLine=-2)

#get rid of some of the choices manually
checkParams<-checkParams[-c(1,2),]
clSmaller<-compareChoices(paramMatrix=checkParams)

## Check subsampling for a couple

#following code takes around 1+ minutes to run because of the subsampling that is redone each time:
system.time(clusterTrack<-compareChoices(simData, ks=5,
                                         alphas=c(0.1), findBestK=c(TRUE,FALSE),
                                         sequential=c(FALSE),
                                         subsample=c(TRUE),
                                         removeSil=c(TRUE,FALSE), clusterMethod=c("pam","tight","hierarchical"),
                                         clusterDArgs = list(minSize = 5,kRange=2:15),ncores=1,random.seed=48120))

#should create error caught by compareChoices.
clusterTrack<-compareChoices(simData, ks=5,
                             alphas=c(0.1), findBestK=c(TRUE,FALSE),
                             sequential=c(FALSE),
                             subsample=c(TRUE),
                             removeSil=c(TRUE,FALSE), clusterMethod=c("pam","tight","hierarchical"),
                             clusterDArgs = list(minSize = 5,kRange=2:15),ncores=1,random.seed=48120))
