# The path of "R-proj"
proj_path<-'/Users/ruibo/Desktop/mvpc/mvpc_v2/MVPC'
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")

source(paste(src_path,'/all_functions.R',sep=""))

logis <- function(x,y){
  return(1/(1+exp(-(x+y))))
}

n = 10000
a =  rnorm(n, 0, 1)
b =  rnorm(n, 0, 1)
e =  rnorm(n, 0, 1)
# 
# a = -n:n
# a = sample(a)
# b = sample(a)


data_t = data.frame(a)
data_t[,2]=b

cor(data_t)
suffStat_t  = list(C=cor(data_t), n=n)
gaussCItest(1, 2, c(), suffStat_t)


r  = a + b >  0


data_t2 = data_t[r,]
cor(data_t2)
suffStat_t2 = list(data=data_t2)
gaussCItest_td(1, 2, c(), suffStat_t2)

c= sum(r)/length(a)

data_t3=data_t2
data_t3[,1] = data_t3[,1]-1000000*exp(-(1/1000000)*(data_t3[,1]+data_t3[,2]))
 
# data_t3[,1]=data_t2[,1]+c*data_t2[,2]

cor(data_t3)
suffStat_t3 = list(C=cor(data_t3), n=length(data_t3[,1]))
gaussCItest(1, 2, c(), suffStat_t3)
# 
par(mfrow=c(1,2))
p_data = data_t2
plot(p_data[,1],p_data[,2], pch=19)
p_data = data_t3
plot(p_data[,1],p_data[,2], pch=19)

