/* Code for doing double for-loop for subsampling
*/

#include <Rcpp.h>
using namespace Rcpp;
#include <sstream>

#include <unordered_map>
#include <string>
#include <iostream>
using namespace std;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
std::string makeHashKey(int i, int j, string delim){
	std::string result;
	if(i<=j)
		result=std::to_string(i)+delim+std::to_string(j);
	else
		result=std::to_string(j)+delim+std::to_string(i);
	return(result);
}

// [[Rcpp::export]]
Rcpp::IntegerVector splitHashKey(string s, string delim){
	Rcpp::IntegerVector token(2);
	int pos = s.find(delim);
	token[0]=std::stoi(s.substr(0, pos));
	s.erase(0, pos + delim.length());
	token[1]=std::stoi(s);

	return token;
}



// [[Rcpp::export]]
Rcpp::IntegerMatrix subsampleLoop(Rcpp::List clList, int N){
	// Create hash of (string of) pairs of indices with count number times subsampled together
	// Will count.
	std::unordered_map<std::string,std::pair< int, int > > countTotal;

	// Preset their bucket size (I *think* this means the number of unique key values, and NOT the size of the values saved in each... need to check.)
	//countTotal.reserve(N*N/2);
	
    Rcpp::IntegerVector clIds;
	Rcpp::IntegerVector clLens;
	Rcpp::List currList;
	std::string delimiter = "-";
	// Interate over the list, i.e. each of the subsamples
	for(int i=0; i< clList.size(); ++i){
		// Should consider whether faster to not input that is list of lists, but rather give two lists, each with vector elements.
		currList=clList[i];
		clIds=currList["clusterIds"];
		clLens=currList["clusterLengths"];
		
		// // Regardless, next part of code would still work, once define those vectors.

		// Calculate cumulative sum.
		int nclusters=clLens.size();
		int nsamples=clIds.size();
		
		Rcpp::IntegerVector cumSumVals=Rcpp::cumsum(clLens);

		// just made for loops to define starts, but should be Rcpp functions to do it
		// Basically want: (0,cumSumVals[Rcpp::Range(0, nclusters-2)])
		Rcpp::IntegerVector indexClusterStarts(nclusters);
		indexClusterStarts[0]=0;
		for(int ii=1; ii<nclusters;++ii){
			indexClusterStarts[ii]=cumSumVals[ii-1];
		}
		
		Rcpp::IntegerVector indexClusterEnds(nclusters);
		for(int ii=0; ii< nclusters-1;++ii){
			indexClusterEnds[ii]=indexClusterStarts[ii+1]-1;
		}
		indexClusterEnds[nclusters-1]=nsamples-1;
		
		// Iterate over each position in the clusterIds vector
		// Add to hash
		for(int j=0; j<nsamples ; j++){
			for(int k=0; k<j;k++){
				std::string hashKey=makeHashKey(clIds[j],clIds[k],delimiter);
				// 1) increase countSampled
				countTotal[hashKey].first+= 1;
				// 2) add only those clustered with j to countClustered[j]
				LogicalVector clusterj = j <= indexClusterEnds & j >= indexClusterStarts;
				LogicalVector clusterk = k <= indexClusterEnds & k >= indexClusterStarts;
				if(is_true(all(clusterj == clusterk))){
					countTotal[hashKey].second+= 1;
				}
			}
		}
	}

	// Add in pairs missing. Iterate over indices (R, so 1-based)
	// For those missing, given them a NA for the total
	for(int i=1; i< N+1;++i){
		for(int j=1; j< i;++j){
			std::string hashKey=makeHashKey(i,j,delimiter);
			if (countTotal.count(hashKey)==0){
				countTotal[hashKey].first=0;
				countTotal[hashKey].second=NA_INTEGER;
			}
			
		}
	}
	
	// Pull out the results
	int nbuckets=countTotal.bucket_count();
	StringVector pairName(nbuckets);
	IntegerVector denom(nbuckets);
	IntegerVector num(nbuckets);
	/* Note: Problem, nbuckets greater than the number of entries because reserved in beginning to try improve memory. 
	Creates too many. 'ind' counts how many real buckets there are*/
	int ind=0;

	IntegerMatrix results(nbuckets,4);
    for ( auto it = countTotal.begin(); it != countTotal.end(); ++it ){
 		if((it->second).first>0 | (it->second).second==NA_INTEGER){
			IntegerVector currPairIndices=splitHashKey(it->first, delimiter);
			results.row(ind)=IntegerVector::create(currPairIndices[0], currPairIndices[1],(it->second).second,(it->second).first);
			
		}
		ind++;
    }
	results=results(Range(0,ind-1),_);
	

	return(results);
}

/* Don't know how to use iterator for IntegerMatrix so that this would work.*/
// // [[Rcpp::export]]
// bool wayToSort(IntegerVector i, IntegerVector j) {
// 	bool out;
// 	if(i[1] == j[1]) out=i[2] < j[2];
// 	else out=i[1]<j[1];
// 	return(out)
// }

/*** Summary of output created from subsampling (input to this function):
	returns a list "DList", each element:
            #each element is list: 
			# 1) "clusterIds": vector of length na.omit(classX) of the original indices, where ids in clusters are adjacent in the vector and 
			# 2) "clusterLengths": another vector of length K indicating length of each cluster (allows to decode where the cluster stopes in the above vector), 
			# What does this do with NAs? Removes them -- not included.
*/

/*** R code for testing
setwd("~/Documents/Sequencing/SingleCell/clusterExperiment")
library(microbenchmark)
library(Rcpp)
#To run default
source("R/subsampleLoop.R")
sourceCpp("src/subsampleLoop.cpp")

#create some sample data
sub1<-list(clusterIds=c(10,2,8,5,1,12),clusterLengths=c(2,3,1))
sub2<-list(clusterIds=c(10,2,3,4,8,6,5,1),clusterLengths=c(3,3,2))
sub3<-list(clusterIds=c(10,3,4,5,1),clusterLengths=c(2,1,2))
exOut<-list(sub1,sub2,sub3)

N<-12

test<-subsampleLoop(exOut, N)
test[order(test[,1],test[,2]),]
Cfun<-function(clusterList,N){
   matResults<-subsampleLoop(clusterList, N)
   ord<-order(matResults[,2],matResults[,1])
	matResults[ord,3]/matResults[ord,4]
}
Rfun<-function(clusterList,N){
	unlist(lapply(2:N,function(jj){searchForPairs(jj,clusterList=clusterList,N=N)}))
}

N<-12
DbarC<-matrix(0,N,N)
DbarR<-matrix(0,N,N)
DbarR[upper.tri(DbarR, diag = FALSE)]<-Rfun(exOut,N)
DbarC[upper.tri(DbarC, diag = FALSE)]<-Cfun(exOut,N)


microbenchmark(
  Rfun(exOut,N),
  Cfun(exOut,N)
)

### Create bigger data
Nbig<-1000
subN<-100
size<-floor(0.7*Nbig)
set.seed(2128950)
clustLens<-sample(10:120,size=10)
clustLens<-c(clustLens,size-sum(clustLens))
if(sum(clustLens)!=size) stop("error in setup")
createFun<-function(){
	ids<-sample(1:Nbig,size,replace=FALSE)
	lens<-sample(clustLens)
	return(list("clusterIds"=ids,clusterLengths=lens))
}
testSubList<-replicate(subN, expr=createFun(), simplify = FALSE)


system.time(testR<-Rfun(testSubList,Nbig))
system.time(testC<-Cfun(testSubList,Nbig))
> system.time(Rfun(testSubList,N))
   user  system elapsed 
 58.702   1.989  61.559 
> system.time(Cfun(testSubList,N))
   user  system elapsed 
 29.070   0.386  29.327 


*/



