args<-commandArgs(TRUE)
if(length(args)==0) stop("Usage should be 'RScript readMemoryLog.R <fileName>' where <fileName> is the name of the file with the memory log.")
file<-args[1]
#file<-"memoryLogger_StandardSmall.txt"
log<-read.table(file,header=FALSE,fill=TRUE,stringsAsFactors=FALSE)
nms<-head(unlist(log[1,]),-1)
log<-log[log[,1]=="Mem:",-1]
names(log)<-nms

#from:
#https://stackoverflow.com/questions/10910688/converting-kilobytes-megabytes-etc-to-bytes-in-r
convb <- function(x){
  ptn <- "(\\d*(.\\d+)*)(.*)"
  num  <- as.numeric(sub(ptn, "\\1", x))
  unit <- sub(ptn, "\\3", x)             
  unit[unit==""] <- "1" 

  mult <- c("1"=1, "K"=1024, "M"=1024^2, "G"=1024^3)
  num * unname(mult[unit])
}
maxMemUsed<-log$used[which.max(convb(log$used))]
cat("Maximum Memory Used:",maxMemUsed,"\n",file=stdout())