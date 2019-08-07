# The path of "R-proj"
proj_path<-'/Users/ruibo/Desktop/mvpc/mvpc_v2/MVPC'
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")
source(paste(src_path,'/all_functions.R',sep=""))

# Data generation
n = 10000
a =  rnorm(n, 0, 1)
b =  rnorm(n, 0, 1)
e =  rnorm(n, 0, 1)
data_t = data.frame(a)
data_t[,2]=b

# Independence test 
suffStat_t  = list(C=cor(data_t), n=n)
gaussCItest(1, 2, c(), suffStat_t)

# Missing values generation 
w <- a + b + e
r  = w >  1

# Independence test 
data_t2 = data_t[r,]
suffStat_t2 = list(data=data_t2)
gaussCItest_td(1, 2, c(), suffStat_t2)

# logistic-reweighted method
rx = as.integer(r) 
logidata <- data.frame(w, rx)
glm.fit <- glm(rx ~ w, data = logidata, family = binomial)
glm.probs = predict(glm.fit, newdata = data.frame(w), type = "response")
beta = sum(rx)/length(rx) * (1/glm.probs)
beta = beta[r]

suffStat_t  = list(C=wtd.cors(data_t2,data_t2,beta), n=sum(rx))
gaussCItest(1, 2, c(), suffStat_t)

# what if we involve unnecessary correction ?

