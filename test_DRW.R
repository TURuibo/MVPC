library(DescTools)
library(mipfp)
## The path of "R-proj"
proj_path<-getwd()
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")
source(paste(src_path,'/all_functions.R',sep=""))

n <- 1000000
# set.seed(123)

x_ <- rbinom(n, 1, pr=1/2)
pz_ <- plogis( 2*x_-1); z_ <- rbinom(n, 1, pr=pz_)
py_ <- plogis( 2*z_-1); y_ <- rbinom(n, 1, pr=py_)
dat <- cbind(x_,y_,z_)

binCItest(1,2,c(3), list(dm = dat, adaptDF = FALSE)) # 0.36, not signif.

# Add missingness mechanism
pw <- plogis( x_ + y_); w_ <- rbinom(n, 1, prob = pw) # {0,1}
pw2 <- plogis( 2*x_ + y_); w2_ <- rbinom(n, 1, prob = pw2) # {0,1}

pr1 <- plogis( w_-0.5); xr_ <- rbinom(n, 1, prob = pr1)
pr <- plogis( 2*w2_-1); yr_ <- rbinom(n, 1, prob = pr)

x = x_
y = y_
z = z_
w = w_
w2 = w2_

xr= xr_
x[xr_==1] = NA
y[yr_==1] = NA

dat = cbind(x,y,z,w,w2)
prt = list()
prt[[1]] = c(4)
prt[[2]] = c(5)

prt_m<-data.frame(m=c(1,2))
prt_m[['prt']]<-prt

suffStat = list(data = dat, prt_m=prt_m, adaptDF = FALSE)

pr = sum(xr_ == 0 & yr_ == 0)/n
pw00 = (sum(xr_ == 0 & w ==0)/sum(w==0))*(sum(yr_ == 0 & w2 ==0)/sum(w2==0))
pw11 = (sum(xr_ == 0 & w ==1)/sum(w==1))*(sum(yr_ == 0 & w2 ==1)/sum(w2==1))
pw01 = (sum(xr_ == 0 & w ==0)/sum(w==0))*(sum(yr_ == 0 & w2 ==1)/sum(w2==1))
pw10 = (sum(xr_ == 0 & w ==1)/sum(w==1))*(sum(yr_ == 0 & w2 ==0)/sum(w2==0))

p=w

p[w== 0 & w2 ==0 ]= pr/pw00
p[w== 0 & w2 ==1 ]= pr/pw01
p[w== 1 & w2 ==0 ]= pr/pw10
p[w== 1 & w2 ==1 ]= pr/pw11

binCItest_td(1, 2, c(3), suffStat)
binCItest_drw(1, 2, c(3), suffStat)




# Test 2 

n <- 1000000
# set.seed(123)

x_ <- rbinom(n, 1, pr=1/2)
pz_ <- plogis( 2*x_-1); z_ <- rbinom(n, 1, pr=pz_)
py_ <- plogis( 2*z_-1); y_ <- rbinom(n, 1, pr=py_)
dat <- cbind(x_,y_,z_)

binCItest(1,2,c(3), list(dm = dat, adaptDF = FALSE)) # 0.36, not signif.

# Add missingness mechanism
pw <- plogis( x_ + y_); w_ <- rbinom(n, 1, prob = pw) # {0,1}
pw2 <- plogis( 2*x_ + y_); w2_ <- rbinom(n, 1, prob = pw2) # {0,1}
pr <- plogis( 2*w2_-1); yr_ <- rbinom(n, 1, prob = pr)

# pr1 <- plogis( w_-0.5); xr_ <- rbinom(n, 1, prob = pr1)

x = x_
y = y_
z = z_
w = w_
w2 = w2_

# x[xr_==1] = NA
y[yr_==1] = NA

dat = cbind(x,y,z,w,w2)

prt_m<-data.frame(m=c(2))
prt = list()
prt[[1]] = c(5)
prt_m[['prt']]<-prt

suffStat = list(data = dat, prt_m=prt_m, adaptDF = FALSE)

binCItest_td(1, 2, c(3), suffStat)
binCItest_drw(1, 2, c(3), suffStat)




# Test 3 

n <- 1000000
set.seed(123)

x_ <- rbinom(n, 1, pr=1/2)
pz_ <- plogis( 2*x_-1); z_ <- rbinom(n, 1, pr=pz_)
py_ <- plogis( 2*z_-1); y_ <- rbinom(n, 1, pr=py_)
dat <- cbind(x_,y_,z_)

binCItest(1,2,c(3), list(dm = dat, adaptDF = FALSE)) # 0.36, not signif.

# Add missingness mechanism
pw <- plogis( x_ + y_); w_ <- rbinom(n, 1, prob = pw) # {0,1}
pr <- plogis( 2*w_-1); yr_ <- rbinom(n, 1, prob = pr)
pr1 <- plogis( w_-0.5); xr_ <- rbinom(n, 1, prob = pr1)

x = x_
y = y_
z = z_
w = w_

x[xr_==1] = NA
y[yr_==1] = NA

dat = cbind(x,y,z,w,w2)

prt_m<-data.frame(m=c(1, 2))
prt = list()
prt[[1]] = c(4)
prt[[2]] = c(4)
prt_m[['prt']]<-prt

suffStat = list(data = dat, prt_m=prt_m, adaptDF = FALSE)

binCItest_td(1, 2, c(3), suffStat)
binCItest_drw(1, 2, c(3), suffStat)
