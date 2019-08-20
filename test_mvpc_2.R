# The path of "R-proj"
proj_path<-getwd()
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")

source(paste(src_path,'/all_functions.R',sep=""))

exp = "mar_cor"
times  = 10
n = 10000
p=20
data_path<-"/Users/ruibo/Desktop/mvpc/mvpc-xyz/causality/R-proj/data"

i=10
load(paste(data_path,"/syn/",exp,n,"_",i,".rda",sep=""))

data_mcar = data_all$data_mcar
data_c = data_all$data_c
suffStat = list(data = data_all$suffStat$data)
myCPDAG=data_all$cpdag
myDAG=data_all$dag

sk<-mvpc(suffStat, gaussCItest_td,PermCCItest, alpha=0.01, p=p)

sk$skel_pre@pMax == sk$skel@pMax

sk$skel_pre@pMax < sk$skel_vs@pMax

res_pc<-pc(suffStat, gaussCItest_td, alpha=0.01, p=p)
# 
# suffStat  = list(C = cor(data_mcar),n=length(data_mcar[,1]))
# res_ref<-pc(suffStat, gaussCItest, alpha=0.01, p=p)
# 
# shd(res_ref, myCPDAG)
# shd(res_mvpc, myCPDAG)
# shd(res_pc, myCPDAG)  

a= data_all$suffStat$sup_var

prt = list()
m = c(1)
prt[[1]]=c(5)


for(i in 1:length(a$m)){
  m[i] = a$m[[i]]
  prt[[i]] = c(a$sup[i])
  print(prt[[i]])
  
}


prt_m<-data.frame(m=m)
prt_m[['prt']]<-prt

data_t2 =  data_all$suffStat$data

suffStat = list(data = data_t2)#, prt_m = prt_m)

PermCCItest(1, 15, c(7), suffStat)
gaussCItest_td(1, 15, c(7), suffStat)


