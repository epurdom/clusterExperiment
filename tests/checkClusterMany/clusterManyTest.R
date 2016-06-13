#Usage: nohup RScript clusterManyTest.R <tagString> <compareTo(optional)> &


library(clusterExperiment) 
load("L5_sumExp.rda")
ncores<-5
args<-commandArgs(TRUE)
if(length(args)==0) stop("Usage should be 'RScript clusterManyTest.R <tagString>' where <tagString> will be name on saved file of output.")
tag<-args[1]
fixedVersion<-if(length(args)==2) args[2] else "fixedClusterManyResult.txt"

x<-sessionInfo()
version<-x$otherPkgs[["clusterExperiment"]][["Version"]]
nm<-paste(tag,"_",version,sep="")

outfile<-paste(nm,".Rout",sep="")
cat("Results for test of",version,"\n",file=outfile)
cat("-------------------\n",file=outfile,append=TRUE)
cat("Running clusterMany...",file=outfile,append=TRUE)
cl <-clusterMany(l5, dimReduce = "PCA", nPCADims = 50, isCount=TRUE,
                 ks=4:8, clusterFunction="hierarchical01",
                 alphas=c(0.2,0.3), subsample=TRUE, sequential=TRUE,
                 ncores=ncores, subsampleArgs=list(resamp.num=20,
                                              clusterFunction="kmeans",
                                              clusterArgs=list(nstart=1)),
                 seqArgs=list(beta=0.9,k.min=3,verbose=FALSE),
                 clusterDArgs=list(minSize=5, verbose=FALSE),
                 random.seed=21321, run=TRUE)
#save(cl, file=paste(tag,"_",version,".rda",sep=""))
cat("done.",file=outfile,append=TRUE)
mat<-clusterMatrix(cl)
row.names(mat)<-colnames(cl)
matFile<-paste(nm,".txt",sep="")
write.table(mat,file=matFile,sep=",",col.names = TRUE,row.names = TRUE)

cat("Current Version:",version,"\n",file=outfile,append=TRUE)
cat("User-given tag:",tag,"\n",file=outfile,append=TRUE)
##Read both in, just to make sure not catching differences due write/read differences
cat("Compare",matFile,"to fixed version (", fixedVersion,")", ":\n",file=outfile,append=TRUE)
compMat<-read.table(fixedVersion,sep=",",header=TRUE)
newMat<-read.table(matFile,sep=",",header=TRUE)
compResult<-all.equal(compMat,newMat)
printResult<-if(isTRUE(compResult)) "Yes" else "No"
cat("Are all entries the same?\n",printResult,"\n",file=outfile,append=TRUE)
cat("-------------------\n",file=outfile,append=TRUE)
cat("Complete Session Info:\n",file=outfile,append=TRUE)
cat(paste(capture.output(x),collapse="\n"),file=outfile,append=TRUE )


