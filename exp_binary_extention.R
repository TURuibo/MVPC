library(DescTools)
## The path of "R-proj"
proj_path<-'/Users/ruibo/Desktop/mvpc/mvpc_v2/MVPC'
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")
source(paste(src_path,'/all_functions.R',sep=""))


## Simulate *independent data of {0,1}-variables:
n <- 10000000
set.seed(123)

x <- rbinom(n, 1, pr=1/2)
y <- rbinom(n, 1, pr=1/2)
# z <- rbinom(n, 1, pr=1/2)
# dat <- cbind(x,y,z)
table1 <- ftable(xtabs(~ x+y))
dat <- cbind(x,y)

binCItest(1,2,c(), list(dm = dat, adaptDF = FALSE)) # 0.36, not signif.
binCItest(1,2,c(), list(dm = dat, adaptDF = TRUE )) # the same, here

# Add missingness mechanism
pw <- plogis( x + y); w <- rbinom(n, 1, prob = pw) # {0,1}
ind <- 1:n
pr <- plogis( 2*w-1); xr <- rbinom(n, 1, prob = pr)
r<-ind[xr==1]
x[r] = NA
w_ = w

x = x[ind[xr==0]]
y = y[ind[xr==0]]
w = w[ind[xr==0]]
xr= xr[ind[xr == 0]]
dat <- cbind(x,y)

table2 <- ftable(xtabs(~ x+y))

## Test-wise deletion doesn't work 
binCItest_td(1,2,c(), list(dm = dat, adaptDF = FALSE))



##  Permutation-based method does not work
ans <- ftable(xtabs(~ x+y+w+xr))
dat <- cbind(x,y,w,xr)

# step 1: get generative model
px1w0 = sum(ans[5]+ans[7])/sum(ans[1]+ ans[3]+ans[5]+ans[7])
px1w1 = sum(ans[6]+ans[8])/sum(ans[2]+ ans[4]+ans[6]+ans[8])
py1w0 = sum(ans[3]+ans[7])/sum(ans[1]+ ans[3]+ans[5]+ans[7])
py1w1 = sum(ans[4]+ans[8])/sum(ans[2]+ ans[4]+ans[6]+ans[8])
  
# step 2: input shuffle w 
w = sample(w_)[1:length(x)]
p_xc = rep(px1w0, length(x))
p_xc[w==1] = px1w1
p_yc = rep(py1w0, length(x))
p_yc[w==1] = py1w1
x_c <- rbinom(length(x), 1, prob = p_xc) # {0,1}
y_c <- rbinom(length(x), 1, prob = p_xc) # {0,1}
dat <- cbind(x_c, y_c)

## Test-wise deletion might introduce more extra edges in PC shown:
binCItest_td(1,2,c(), list(dm = dat, adaptDF = FALSE))
table3 <-  ftable(xtabs(~ x_c+y_c))

GTest(a,correct="none")            # "none" "williams" "yates"


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
