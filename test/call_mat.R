# The path of "R-proj"
proj_path<-getwd()
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")

source(paste(src_path,'/all_functions.R',sep=""))

#************ set up the environment of R.Matlab ************

library(R.matlab)
options(matlab="/Applications/MATLAB_R2019b.app/bin/matlab")
Matlab$startServer()

#************ set up the environment of R.Matlab ************

matlab <- Matlab()
isOpen <- open(matlab)
path = "/Users/ruibo/Desktop/mvpc/mvpc_v2/MVPC/src/KCI-test"
addpath = paste0("addpath('",path,"')", sep = "")
evaluate(matlab, addpath)
#************ END: set up the environment of R.Matlab ************

x <- rnorm(1000, mean=0, sd=1)
y <- rnorm(1000, mean=0, sd=1)
w <- x + y + rnorm(1000, mean=0, sd=1)

data = cbind(x,y)
suffStat=list(C= cor(data),n=length(data[,1]))
gaussCItest(1, 2, c(), suffStat)

ind = 1:length(data[,1])
data1= data[ind[w>-1], ]
suffStat1=list(C= cor(data1),n=length(data1[,1]))
gaussCItest(1, 2, c(), suffStat1)

## only support one parent in the current version
beta = k.weights(w[w>-1],w)
suffStat2=list(C= wtd.cors(data1,data1,beta),n=length(beta))
gaussCItest(1, 2, c(), suffStat2)

close(matlab)
