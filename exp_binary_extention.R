library(DescTools)
library(mipfp)
## The path of "R-proj"
proj_path<-getwd()
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")
source(paste(src_path,'/all_functions.R',sep=""))

## Data generation
n <- 1000000
set.seed(123)

x_ <- rbinom(n, 1, pr=1/2)
pz_ <- plogis( 2*x_-1); z_ <- rbinom(n, 1, pr=pz_)
py_ <- plogis( 2*z_-1); y_ <- rbinom(n, 1, pr=py_)
dat <- cbind(x_,y_,z_)

binCItest(1,2,c(3), list(dm = dat, adaptDF = FALSE)) # 0.36, not signif.

# Add missingness mechanism
pw <- plogis( x_ + y_); w_ <- rbinom(n, 1, prob = pw) # {0,1}
pr <- plogis( 2*w_-1); xr_ <- rbinom(n, 1, prob = pr)

x = x_[xr_==0]
y = y_[xr_==0]
z = z_[xr_==0]
w = w_[xr_==0]
xr= xr_[xr_==0]

## Test-wise deletion doesn't work 
dat = cbind(x,y,z)
binCItest_td(1,2,c(3), list(data = dat, adaptDF = FALSE))

## Permutation-based Correction Method works
dat <- cbind(x, y, z)

dat_w0 <- dat[w==0,]
dat_w1 <- dat[w==1,]

p.xy.w0 <- apply(dat_w0,2, function(x) sum(x)/length(x))
cor.xy.w0 <- cor(dat_w0)
p.joint.w0 <- ObtainMultBinaryDist(corr = cor.xy.w0, marg.probs = p.xy.w0 )

p.xy.w1 <- apply(dat_w1 ,2, function(x) sum(x)/length(x))
cor.xy.w1 <- cor(dat_w1 )
p.joint.w1 <- ObtainMultBinaryDist(corr = cor.xy.w1, marg.probs = p.xy.w1)

w = sample(w_)[1:length(x)]
data_w0 <- RMultBinary(n = sum(w==0), mult.bin.dist = p.joint.w0)$binary.sequences
data_w1 <- RMultBinary(n = sum(w==1), mult.bin.dist = p.joint.w1)$binary.sequences
dat <- rbind(data_w0,data_w1 )
binCItest(1,2,c(3), list(dm = dat , adaptDF = FALSE))

## Permutation-based Correction Test works
x = x_
y = y_
z = z_
w = w_

x[xr_==1] = NA
data = cbind(x,y,z,w)

prt = list()
m = c(1)
prt[[1]]=c(4)

prt_m<-data.frame(m=m)
prt_m[['prt']]<-prt

suffStat = list(data = data , prt_m = prt_m)
binCItest_td(1,2,c(3), suffStat)
binPermCCItest(1,2,c(3), suffStat)

## Successful: DRW correction method works for two variables case 
pr0 = sum(xr == 0)/ length(xr)
pr0w0 = sum(xr == 0 & w==0)/sum( w==0 )
pr0w1 = sum(xr == 0 & w==1)/sum( w==1 )
c0 = pr0/pr0w0
c1 = pr0/pr0w1
x = x[xr==0]
y = y[xr==0]
w = w[xr==0]
Data.xtabs0 = c0 * xtabs( ~ x[w==0]+ y[w==0])
Data.xtabs1 = c1 * xtabs( ~ x[w==1]+ y[w==1])
Data.xtabs = Data.xtabs0  + Data.xtabs1
GTest(Data.xtabs,correct="none")            # "none" "williams" "yates"

##  Permutation-based method does not work
# step 1: get generative model
px1w0 = sum((dat[,'w']==0 )&(dat[,'x']==1 ))/sum(dat[,'w']==0 )
px1w1 = sum((dat[,'w']==1 )&(dat[,'x']==1 ))/sum(dat[,'w']==1 )
py1w0 = sum((dat[,'w']==0 )&(dat[,'y']==1 ))/sum(dat[,'w']==0 )
py1w1 = sum((dat[,'w']==1 )&(dat[,'y']==1 ))/sum(dat[,'w']==1 )

# step 2: input shuffle w 
w = sample(w_)[1:length(x)]
p_xc = rep(px1w0, length(x))
p_xc[w==1] = px1w1

p_yc = rep(py1w0, length(x))
p_yc[w==1] = py1w1

x_c <- rbinom(length(x), 1, prob = p_xc) # {0,1}
y_c <- rbinom(length(x), 1, prob = p_yc) # {0,1}

dat <- cbind(x_c, y_c,w)
binCItest(1,2,c(), list(dm = dat , adaptDF = FALSE))
