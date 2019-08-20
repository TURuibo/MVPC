# The path of "R-proj"
proj_path<-getwd()
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")

source(paste(src_path,'/all_functions.R',sep=""))

exp_10000 <- comp_com_td_mvpc(n=10000, "mnar_cor")
exp_1000 <-comp_com_td_mvpc(n=1000,"mnar_cor")
exp_100 <-comp_com_td_mvpc(n=100,"mnar_cor")

result<-list(exp_100,exp_1000,exp_10000)

# Mean and standard deviation
mn<-c()
std<-c()
for(s in result){
  mn1<-array(unlist(lapply(s[1:3],mean)))
  std1<-array(unlist(lapply(s[1:3],sd)))
  mn<-c(mn, mn1)
  std<-c(std, std1)
}

