ncores<-5
large<-TRUE
resamp.num<-200
tag<-"memTest_LargeBig"

#Usage: nohup RScript clusterManyTest.R <tagString> <compareTo(optional)> &
# If get that corrupted file, probably copied from laptop or elsewhere that only has tag
# Do git lfs checkout L5_sumExp.rda
library(devtools)
library(profmem)
load_all()
load("L5_sumExp.rda")
outpath<-"memTestDirectory"
if(!file.exists(outpath)) dir.create(outpath)
fixedVersion<-"fixedClusterManyResult.txt"

x<-sessionInfo()
version<-x$otherPkgs[["clusterExperiment"]][["Version"]]
nm<-paste(tag,"_",version,sep="")

outfile<-file.path(outpath,paste(nm,".Rout",sep=""))
cat("Results for test of",version,"\n",file=outfile)
cat("-------------------\n",file=outfile,append=TRUE)
cat("Running clusterMany...",file=outfile,append=TRUE)
# Old version: 
# cl <-clusterMany(l5, reduceMethod = "PCA", nReducedDims = 50, isCount=TRUE,
#                  ks=4:8, clusterFunction="hierarchical01",
#                  alphas=c(0.2,0.3), subsample=TRUE, sequential=TRUE,
#                  ncores=ncores, subsampleArgs=list(resamp.num=20,
#                                               clusterFunction="kmeans",
#                                               clusterArgs=list(nstart=1)),
#                  seqArgs=list(beta=0.9,k.min=3,verbose=FALSE),
#                  mainClusterArgs=list(minSize=5, verbose=FALSE),
#                  random.seed=21321, run=TRUE)
sttm<-proc.time()
cl <-clusterMany(l5, reduceMethod = "PCA", nReducedDims = 50, isCount=TRUE,
                 ks=4:8, clusterFunction="hierarchical01",
                 beta=0.9, minSize=5, mainClusterArgs=list(clusterArgs=list("whichHierDist"="dist")), #added this to be back-compatible with previous defauls.
				 seqArgs=list(top.can=15),#added this to be back-compatible with previous defauls.
                 alphas=c(0.2,0.3), subsample=TRUE, sequential=TRUE,
                 ncores=ncores, subsampleArgs=list(resamp.num=resamp.num,largeDataset=large,
                                                   clusterFunction="kmeans",
                                                   clusterArgs=list(nstart=1)),
                 random.seed=21321, run=TRUE)

endtm<-proc.time()
tm<-endtm-sttm
#save(cl, file=paste(tag,"_",version,".rda",sep=""))
cat("done.\n",file=outfile,append=TRUE)
cat(paste("Ellapsed Time:",tm[3]/60,"minutes\n"),file=outfile,append=TRUE)
print(tm) #just in case.
