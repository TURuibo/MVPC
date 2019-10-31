proj_path<-getwd()
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")

source(paste(src_path,'/all_functions.R',sep=""))

## ********* Synthethic Binary Data Generation ********* 
data.name = 'df_sim_mar_0.9_num'
data.file <- paste(data_path,'/',data.name,'.csv',sep="")
MyData <- read.csv(file=data.file, header=TRUE, sep=",")
MyData$X <- NULL 

## ********* Correction ********* 
suffStat = list(data=MyData[1:1000, ],adaptDF=FALSE)
res_mvpc<-mvpc(suffStat, binCItest_td, binPermCCItest, alpha=0.05, p=222)
save(res_mvpc, file = "res_mvpc.RData")

## ********* Test-Wise Deletion ********* 
suffStat = list(data=MyData[1:1000, ],adaptDF=FALSE)
sample_size <<- c()
res_td_pc<-pc(suffStat, binCItest_td, alpha=0.05, p=222)
save(res_td_pc, file = "res_td_pc.RData")
sample_size <- mean(sample_size)

# Complete data evaluation
data.name = 'df_sim_no_mar_0.9_num'
data.file <- paste(data_path,'/',data.name,'.csv',sep="")
MyData <- read.csv(file=data.file, header=TRUE, sep=",")
MyData$X <- NULL 

suffStat = list(dm=MyData[1:1000, ],adaptDF=FALSE)
res_pc<-pc(suffStat, binCItest, alpha=0.05, p=222)
res_pc_full<-res_pc
save(res_pc_full, file = "res_pc_full.RData")

suffStat = list(dm=MyData[1:sample_size, ],adaptDF=FALSE)
res_pc_ref<-pc(suffStat, binCItest, alpha=0.05, p=222)
save(res_pc_ref, file = "res_pc_ref.RData")
