# The path of "R-proj"
proj_path<-'/Users/ruibo/Desktop/mvpc/mvpc_v2/MVPC'
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")
source(paste(src_path,'/all_functions.R',sep=""))

mis_f<-function(y){
  if(y < 0) 
    rbinom(1,1,0.99) == 1
  else 
    rbinom(1,1,0.1) == 1
} 

# Data generation
n = 10000
c= rnorm(n, 0, 1)
a = c+rnorm(n, 0, 1)
b =  c+rnorm(n, 0, 1)
e = rnorm(n, 0, 1)
data_t = data.frame(a)
data_t[,2]=b

par(mfrow=c(2,2))
plot(data_t,xlim = c(-4,4),asp=1)
plot(data_t[a>0,],xlim = c(-4,4),asp=1)
plot(data_t[a-0.3*b>0,],xlim = c(-4,4),asp=1)
plot(data_t[a-0.3*b+e>0,],xlim = c(-4,4),asp=1)


# Independence test 
suffStat_t  = list(C=cor(data_t), n=n)
gaussCItest(1, 2, c(), suffStat_t)

# Missing values generation 
w <- a + b + e
# r  = w  >  1
input = w  + rnorm(n, 0, 1)
r<-sapply(input, function(input) mis_f(input))
data_t[,3]=w

# what if we involve unnecessary correction ?
# # Independence test 
data_t2 = data_t
data_t2 = data_t[r,]
suffStat_t1 = list(data=data_t2)
gaussCItest_td(1, 2, c(), suffStat_t1)



beta = w
beta[w<0] = 1/0.99
beta[!w<0] = 1/0.1
beta = sum(rx)/length(rx) * beta
beta = beta[r]

suffStat_t2  = list(C=wtd.cors(data_t2,data_t2,beta), n=sum(rx))
gaussCItest(1, 2, c(), suffStat_t2)


# density estimation 

fhat <- kde(x=w)
fhat_r <- kde(x=w[r])
beta3<- predict(fhat, x=w[r])/predict(fhat_r, x=w[r])
suffStat_t3  = list(C=wtd.cors(data_t2,data_t2,beta3), n=sum(rx))
gaussCItest(1, 2, c(), suffStat_t3)


# DRW independence test 
data_t3 = data_t
data_t3[!r,1] = NA

m = 1
prt=list()
prt[[1]] = c(3)

prt_m<-data.frame(m=m)
prt_m[['prt']]<-prt
suffStat4 <- list(data= data_t3, prt_m= prt_m)

# DRW independence test 
data_t5 = data_t
data_t5[!r,1] = NA
data_t5[,1] = log(1+exp(w))/(sum(r)/length(r))
suffStat5 <- list(C= cor(data_t5), n=length(data_t5[,1]))
DRWCItest(1, 2, c(), suffStat4)

# logistic-reweighted method
# rx = as.integer(r) 
# logidata <- data.frame(rx,w)
# glm.fit <- glm(rx ~  w, data = logidata, family = binomial)
# glm.probs = predict(glm.fit, newdata = data.frame(w), type = "response")
# beta1 = sum(rx)/length(rx) * (1/glm.probs)
# beta1 = beta1[r]
# 
# # Independence test 
# data_t2 = data_t
# data_t2 = data_t[r,]
# suffStat_t2 = list(data=data_t2)
# gaussCItest_td(1, 2, c(), suffStat_t2)
