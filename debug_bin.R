proj_path<-getwd()
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")

source(paste(src_path,'/all_functions.R',sep=""))

set.seed(777)
## ********* Synthethic Binary Data Generation ********* 
i_g = 1
print(paste('graph=', i_g))
data.name = paste('data',i_g,sep='')
graph.name = paste('graph',i_g,sep='')

data.file <- paste(data_path,'/',data.name,'.txt',sep="")
graph.file<- paste(data_path,'/',graph.name,'.txt',sep="")

data = load_bin_data(data.file)
DAG = load_bin_graph(graph.file)
CPDAG = dag2cpdag(DAG)

# Detect colliders and  Colliders' parents
cldr <- detect_colliders(DAG)
cldr_prt <- detect_colliders_prt(DAG, cldr)

# Choose missingness inidcator and their parents
p_m <- create_m_ind(cldr,cldr_prt)
ms = p_m$ms
prt_ms = p_m$prt_ms

prt_m<-data.frame(m=ms)

prt = list()
for(i in 1:length(prt_ms)){
  prt[[i]] = c(prt_ms[i])  
}

prt_m[['prt']]<-prt
prt_m
# Generate missing values 
mask = data != data

for(i in 1:length(ms)){
  nsample = nrow(data)
  p_1 = runif(1,0.1,0.5)
  m_ind = rbinom(nsample, 1, p_1)==1
  m_ind = (data[,prt_ms[i]] == 1) & m_ind
  mask[,ms[i]] = m_ind
}
data_m = data
data_m[mask] = NA

## ********* Correction ********* 
suffStat = list(data=data_m,adaptDF=FALSE)
res_mvpc<-mvpc(suffStat, binCItest_td, binCItest.drw,prt_m, alpha=0.05, p=20)
shd(res_mvpc,CPDAG)
## ********* Test-Wise Deletion ********* 
sample_size <<- c()
suffStat = list(data=data_m,adaptDF=FALSE)
res_td_pc<-pc(suffStat, binCItest_td, alpha=0.05, p=20)
shd(res_td_pc,CPDAG)

## ********* Complete data evaluation *********
suffStat = list(data=data,adaptDF=FALSE)
res_pc<-pc(suffStat, binCItest_td, alpha=0.05, p=20)
shd(res_pc,CPDAG)
