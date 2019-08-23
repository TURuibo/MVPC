proj_path<-getwd()
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")

source(paste(src_path,'/all_functions.R',sep=""))

## ********* Synthethic Binary Data Generation ********* 
data.name = 'data6'
graph.name = 'graph6'

data.file <- paste(data_path,'/',data.name,'.txt',sep="")
graph.file<- paste(data_path,'/',graph.name,'.txt',sep="")

data = load_bin_data(data.file)
CPDAG = load_bin_graph(graph.file)
# Complete data evaluation
suffStat = list(dm=data,adaptDF=FALSE)
res_pc<-pc(suffStat, binCItest, alpha=0.01, p=20)
shd(res_pc,CPDAG)

## ********* Missingness generation ********* 


## ********* Test-Wise Deletion ********* 

## ********* Correction ********* 
res_mvpc<-mvpc(suffStat, binCItest_td, binPermCCItest, alpha=0.01, p=20)




