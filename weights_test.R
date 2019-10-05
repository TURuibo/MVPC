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

checkp1 = unique(p)

checkp1
sum(compute_weights(1, 2, c(3), suffStat)-p)

compute_weights<- function(x, y, S, suffStat){
  data = suffStat$data
  n.sample = dim(data)[1]
  weights = rep(0, n.sample)
  
  # Detection of parents of missingness indicators
  ind_test <- c(x, y, S)
  ind_W <- get_prt_m_xys(c(x,y,S), suffStat)  # Get parents the {xyS} missingness indicators
  ind_W <- setdiff(ind_W,ind_test)
  if(length(ind_W)==0){return(binCItest_td(x,y,S,suffStat))}
  
  ind_W <- c(ind_W, get_prt_m_xys(ind_W, suffStat) ) # Get parents the W missingness indicators
  ind_W <- setdiff(ind_W,ind_test)
  
  # Get ri and corresponding wi
  rw = get_rw_pair(x,y,S,ind_W,suffStat)
  # Get the weights check table weights <==> W
  num.W <- length(ind_W)
  comb.W <- bincombinations(length(ind_W))
  weights_tab = weights_check_table(rw,suffStat,ind_W,comb.W)
  # apply the weights check table to assigning weights for each data point
  W = matrix(data = data[,ind_W], nrow = n.sample)
  
  for(i in 1: dim(comb.W)[1]){
    ind = rep(TRUE, n.sample)
    for(j in 1:dim(comb.W)[2]){
      ind = (W[,j] == comb.W[i,j]) & ind
    }
    weights[ind] =  weights_tab[i] 
  }
  return(weights)
}

weights_check_table<-function(rw,suffStat,indW, comb.W){
  weights_tab = c()
  pR <- f_R(rw$r,suffStat)
  n.comb.W <- dim(comb.W)[1]
  for(i in 1:n.comb.W){
    weighti = pR
    for(count in 1:length(rw$r)){
      ri = rw$r[count]
      wi =rw$w[[count]]
      pos_wi = (wi == indW) %*% 1:length(indW)
      val_wi = comb.W[i, pos_wi]
      weighti = weighti / f_weights(ri, wi, val_wi, suffStat)
    }
    weights_tab = c(weights_tab, weighti)
  }
  weights_tab
}
