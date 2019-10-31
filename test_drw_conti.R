#************ set up the environment of R.Matlab ************
library(R.matlab)
close(matlab)
print("Setting up matlab")
options(matlab="/Applications/MATLAB_R2019b.app/bin/matlab")
Matlab$startServer()
matlab <- Matlab()
isOpen <- open(matlab)
path = "/Users/ruibo/Desktop/mvpc/mvpc_v2/MVPC/src/KCI-test"
addpath = paste0("addpath('",path,"')", sep = "")
evaluate(matlab, addpath)
#************ END: set up the environment of R.Matlab ************
## The path of "R-proj"
proj_path<-getwd()
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")
source(paste(src_path,'/all_functions.R',sep=""))


n <- 1000
# set.seed(123)
x <- rnorm(n, mean=0, sd=1)
z <- rnorm(n, mean=0, sd=1) + x
y <- rnorm(n, mean=0, sd=1) + z
w <- x + y + rnorm(n, mean=0, sd=1)
w2 <- 2*x + y + rnorm(n, mean=0, sd=1)
dat <- cbind(x,y,z,w,w2)

suffStat=list(C= cor(dat),n=length(dat[,1]))
gaussCItest(1, 2, c(3), suffStat)

# Add missingness mechanism
dat.m = dat
dat.m[w>-1 , 1] = NA # W -- > Rx

prt = list()
prt[[1]] = c(4)
# prt[[2]] = c(5)
prt_m<-data.frame(m=c(1))
prt_m[['prt']]<-prt

suffStat = list(data = dat.m, prt_m=prt_m, adaptDF = FALSE)
gaussCItest_td(1, 2, c(3), suffStat)
gaussCItest.drw(1, 2, c(3), suffStat)
# gaussCItest.drw(1, 2, c(), suffStat)


# **************** Test 2 **************** 

# Add missingness mechanism
dat.m2 = dat
dat.m2[w>-1 , 1] = NA # W -- > Rx
dat.m2[w2>-1 , 2] = NA # W2 -- > Ry

prt = list()
prt[[1]] = c(4)
prt[[2]] = c(5)
prt_m<-data.frame(m=c(1,2))
prt_m[['prt']]<-prt

suffStat = list(data = dat.m2, prt_m=prt_m, adaptDF = FALSE)
gaussCItest_td(1, 2, c(3), suffStat)
gaussCItest.drw(1, 2, c(3), suffStat)
# **************** END: Test 2 ****************

# **************** Test 3 ****************
# Add missingness mechanism
dat.m3 = dat
dat.m3[w>-1 , 1] = NA # W -- > Rx
dat.m3[w>-1 , 2] = NA # W2 -- > Ry

prt = list()
prt[[1]] = c(4)
prt[[2]] = c(4)
prt_m<-data.frame(m=c(1,2))
prt_m[['prt']]<-prt

suffStat = list(data = dat.m3, prt_m=prt_m, adaptDF = FALSE)
gaussCItest_td(1, 2, c(3), suffStat)
gaussCItest.drw(1, 2, c(3), suffStat)
# **************** END: Test 3 ****************
close(matlab)
