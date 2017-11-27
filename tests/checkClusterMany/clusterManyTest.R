#Usage: nohup RScript clusterManyTest.R <tagString> <compareTo(optional)> &
# If get that corrupted file, probably copied from laptop or elsewhere that only has tag
# Do git lfs checkout L5_sumExp.rda
library(devtools)
library(profmem)
load_all()
#install.packages(pkgs="../../../clusterExperiment",repos=NULL,type="source")
#library(clusterExperiment)
load("L5_sumExp.rda")
outpath<-"resultsDirectory"
if(!file.exists(outpath)) dir.create(outpath)
ncores<-5
args<-commandArgs(TRUE)
if(length(args)==0) stop("Usage should be 'RScript clusterManyTest.R <tagString>' where <tagString> will be name on saved file of output.")
tag<-args[1]
fixedVersion<-if(length(args)==2) args[2] else "fixedClusterManyResult.txt"

x<-sessionInfo()
version<-x$otherPkgs[["clusterExperiment"]][["Version"]]
nm<-paste(tag,"_",version,sep="")
# library(benchmarkme) #to get information about RAM and CPUs on machine:
# print(get_ram())
# print(get_cpu())
# print(sessionInfo())

outfile<-file.path(outpath,paste(nm,".Rout",sep=""))
cat("Results for test of",version,"\n",file=outfile)
cat("-------------------\n",file=outfile,append=TRUE)
cat("Running clusterMany...",file=outfile,append=TRUE)
# Old version: 
# cl <-clusterMany(l5, dimReduce = "PCA", nPCADims = 50, isCount=TRUE,
#                  ks=4:8, clusterFunction="hierarchical01",
#                  alphas=c(0.2,0.3), subsample=TRUE, sequential=TRUE,
#                  ncores=ncores, subsampleArgs=list(resamp.num=20,
#                                               clusterFunction="kmeans",
#                                               clusterArgs=list(nstart=1)),
#                  seqArgs=list(beta=0.9,k.min=3,verbose=FALSE),
#                  mainClusterArgs=list(minSize=5, verbose=FALSE),
#                  random.seed=21321, run=TRUE)
sttm<-proc.time()
cl <-clusterMany(l5, dimReduce = "PCA", nPCADims = 50, isCount=TRUE,
                 ks=4:8, clusterFunction="hierarchical01",
                 beta=0.9, minSize=5, mainClusterArgs=list(clusterArgs=list("whichHierDist"="dist")), #added this to be back-compatible with previous defauls.
				 seqArgs=list(top.can=15),#added this to be back-compatible with previous defauls.
                 alphas=c(0.2,0.3), subsample=TRUE, sequential=TRUE,
                 ncores=ncores, subsampleArgs=list(resamp.num=20,largeDataset=FALSE,
                                                   clusterFunction="kmeans",
                                                   clusterArgs=list(nstart=1)),
                 random.seed=21321, run=TRUE)
endtm<-proc.time()
tm<-endtm-sttm
#save(cl, file=paste(tag,"_",version,".rda",sep=""))
cat("done.\n",file=outfile,append=TRUE)
cat(paste("Ellapsed Time:",tm[3]/60,"minutes\n"),file=outfile,append=TRUE)
cat(paste("Number of clusters:",nClusterings(cl),"\n"),file=outfile,append=TRUE)
cat(paste("Number of genes of Assay:",nrow(cl),"\n"),file=outfile,append=TRUE)
cat(paste("Number of samples of Assay:",ncol(cl),"\n"),file=outfile,append=TRUE)


mat<-clusterMatrix(cl)
row.names(mat)<-colnames(cl)
matFile<-paste(nm,".txt",sep="")
write.table(mat,file=file.path(outpath,matFile),sep=",",col.names = TRUE,row.names = TRUE)

cat("Current Version:",version,"\n",file=outfile,append=TRUE)
cat("User-given tag:",tag,"\n",file=outfile,append=TRUE)
##Read both in, just to make sure not catching differences due write/read differences
cat("Compare",matFile,"to fixed version (", fixedVersion,")", ":\n",file=outfile,append=TRUE)
compMat<-read.table(fixedVersion,sep=",",header=TRUE)
newMat<-read.table(file.path(outpath,matFile),sep=",",header=TRUE)
compResult<-all.equal(compMat,newMat)
printResult<-if(isTRUE(compResult)) "Yes" else "No"
cat("Are all entries the same?\n",printResult,"\n",file=outfile,append=TRUE)
cat("-------------------\n",file=outfile,append=TRUE)
cat("Complete Session Info:\n",file=outfile,append=TRUE)
cat(paste(capture.output(x),collapse="\n"),file=outfile,append=TRUE )


