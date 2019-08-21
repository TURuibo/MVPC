library(DescTools)
library(mipfp)
## The path of "R-proj"
proj_path<-getwd()
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")
source(paste(src_path,'/all_functions.R',sep=""))


## Simulate *independent data of {0,1}-variables:
n <- 1000000
set.seed(123)

x_ <- rbinom(n, 1, pr=1/2)
y_ <- rbinom(n, 1, pr=1/2)
table1 <- ftable(xtabs(~ x+y))

dat <- cbind(x_,y_)

binCItest(1,2,c(), list(dm = dat, adaptDF = FALSE)) # 0.36, not signif.
binCItest(1,2,c(), list(dm = dat, adaptDF = TRUE )) # the same, here

# Add missingness mechanism
pw <- plogis( x_ + y_); w_ <- rbinom(n, 1, prob = pw) # {0,1}
pr <- plogis( 2*w_-1); xr_ <- rbinom(n, 1, prob = pr)
dat_ <- cbind(x_,y_,w_,xr_)
binCItest(1,2,c(4), list(dm = dat_, adaptDF = FALSE)) # 0.36, not signif.
binCItest(1,4,c(3), list(dm = dat_, adaptDF = FALSE)) # 0.36, not signif.


x = x_[xr_==0]
y = y_[xr_==0]
w = w_[xr_==0]
xr= xr_[xr_==0]

dat <- dat_[xr_==0,]

## Test-wise deletion doesn't work 
binCItest_td(1,2,c(), list(dm = dat, adaptDF = FALSE))
binCItest_td(1,2,c(3), list(dm = dat, adaptDF = FALSE))

dat <- cbind(x,y,w,xr)
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

# Why it doesn't work ? because we only care about the marginal distribution here, it has no difference with the x <- w ->  y;
binCItest(1,2,c(3), list(dm = dat, adaptDF = FALSE))

# from Qaqish et al. (2012)
or <- matrix(c(Inf, 0.281, 2.214, 2.214,
               0.281, Inf, 2.214, 2.214,
               2.214, 2.214, Inf, 2.185,
               2.214, 2.214, 2.185, Inf), nrow = 4, ncol = 4, byrow = TRUE)
rownames(or) <- colnames(or) <- c("Parent1", "Parent2", "Sibling1", "Sibling2")

# hypothetical marginal probabilities
p <- c(0.2, 0.4, 0.6, 0.8)

# estimating the joint-distribution
p.joint <- ObtainMultBinaryDist(odds = or, marg.probs = p)

# simulating 100,000 draws from the obtained joint-distribution
y.sim <- RMultBinary(n = 1e5, mult.bin.dist = p.joint)$binary.sequences

# checking results
cat('dim y.sim =', dim(y.sim)[1], 'x', dim(y.sim)[2], '\n')
cat('Estimated marginal probs from simulated data\n')
apply(y.sim,2,mean)
cat('True probabilities\n')
print(p)
cat('Estimated correlation from simulated data\n')
cor(y.sim)
cat('True correlation\n')
Odds2Corr(or,p)$corr

# generating binary outcomes with outcome different than 0, 1
RMultBinary(n = 10, mult.bin.dist = p.joint, 
            target.values = list(c("A", "B"), c(0, 1), c(1, 2), c(100, 101)))



ftable(xtabs(~ x_+y_+w_))
sum((dat_[,'w_']==0 )&(dat_[,'x_']==1 ))/sum(dat_[,'w_']==0 )
sum((dat_[,'w_']==0 )&(dat_[,'x_']==0 ))/sum(dat_[,'w_']==0 )
sum((dat_[,'w_']==0 )&(dat_[,'y_']==1 ))/sum(dat_[,'w_']==0 )
sum((dat_[,'w_']==0 )&(dat_[,'y_']==0 ))/sum(dat_[,'w_']==0 )

sum((dat_[,'w_']==1 )&(dat_[,'y_']==1 ))/sum(dat_[,'w_']==1 )
sum((dat_[,'w_']==1 )&(dat_[,'y_']==0 ))/sum(dat_[,'w_']==1 )
sum((dat_[,'w_']==1 )&(dat_[,'x_']==1 ))/sum(dat_[,'w_']==1 )
sum((dat_[,'w_']==1 )&(dat_[,'x_']==0 ))/sum(dat_[,'w_']==1 )




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
