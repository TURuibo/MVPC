proj_path<-getwd()
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")

source(paste(src_path,'/all_functions.R',sep=""))

set.seed(100)
td_pc = list()
mvpc_drw = list()
mvpc_permc = list()
ref = list()
pc = list()

std_td_pc = list()
std_mvpc_drw = list()
std_mvpc_permc = list()
std_ref = list()
std_pc = list()

r_td_pc = list()
r_mvpc_drw = list()
r_mvpc_permc = list()
r_ref = list()
r_pc = list()

## ********* Synthethic Binary Data Generation ********* 

shd_mvpc_permc = c()
shd_mvpc_drw = c()

shd_ref = c()
shd_pc = c()
shd_td_pc = c()

rp_mvpc_drw = list()
rp_mvpc_permc = list()

rp_ref = list()
rp_pc = list()
rp_td_pc = list()

sz = 600000
set.seed(100)
count = 1
for(i_g in 1:10){
  print(paste('graph=', i_g+10,'number of samples: ', sz))
  data.name = paste('data',i_g+10,sep='')
  graph.name = paste('graph',i_g+10,sep='')
  
  data.file <- paste(data_path,'/',data.name,'.txt',sep="")
  graph.file<- paste(data_path,'/',graph.name,'.txt',sep="")
  
  data = load_bin_data(data.file)
  data = data[1:sz,]
  DAG = load_bin_graph(graph.file)
  CPDAG = dag2cpdag(DAG)
  
  # Detect colliders and  Colliders' parents
  cldr <- detect_colliders(DAG)
  cldr_prt <- detect_colliders_prt(DAG, cldr)
  # Choose missingness inidcator and their parents

  
  p_m <- create_mar_ind(cldr,cldr_prt,
                         num_var=20, 
                         num_extra_e=5, 
                         num_m = 10) 
  ms = p_m$ms
  prt_ms = p_m$prt_ms
  
  # Generate missing values 
  mask = data != data
  
  for(i in 1:length(ms)){
    nsample = nrow(data)
    pr1 <- plogis(2*data[,prt_ms[i]]-0.5); 
    r <- rbinom(nsample, 1, prob = pr1)==1
    mask[,ms[i]] = r
  }
  data_m = data
  data_m[mask] = NA
  prt_m<-data.frame(m=ms)
  
  prt = list()
  for(i in 1:length(prt_ms)){
    prt[[i]] = c(prt_ms[i])  
  }
  
  prt_m[['prt']]<-prt
  ## ********* Correction *********
  suffStat = list(data=data_m,adaptDF=FALSE)
  res_mvpc_permc<-mvpc(suffStat, binCItest_td, binCItest.permc,prt_m, alpha=0.05, p=20)
  
  ## ********* Correction *********
  suffStat = list(data=data_m,adaptDF=FALSE)
  res_mvpc_drw<-mvpc(suffStat, binCItest_td, binCItest.drw, prt_m, alpha=0.05, p=20)
  
  ## ********* Complete data evaluation *********
  suffStat = list(data=data,adaptDF=FALSE)
  res_pc<-pc(suffStat, binCItest_td, alpha=0.05, p=20)
  
  ## ********* Test-Wise Deletion *********
  sample_size <<- c()
  suffStat = list(data=data_m,adaptDF=FALSE)
  res_td_pc<-pc(suffStat, binCItest_td_ref, alpha=0.05, p=20)
  sample_size <- mean(sample_size)
  
  ## ********* MCAR Complete data evaluation *********
  suffStat = list(dm=data[1:sample_size,],adaptDF=FALSE)
  res_mcar_pc<-pc(suffStat, binCItest, alpha=0.05, p=20)
  
  i_ind = i_g
  shd_mvpc_permc[i_ind] =shd(res_mvpc_permc,CPDAG)
  shd_mvpc_drw[i_ind] =shd(res_mvpc_drw,CPDAG)
  shd_ref[i_ind] =shd(res_mcar_pc,CPDAG)
  shd_pc[i_ind] = shd(res_pc,CPDAG)
  shd_td_pc[i_ind] =shd(res_td_pc,CPDAG)
  
  rp_mvpc_permc[[i_ind]]= test_adj(res_mvpc_permc, CPDAG)
  rp_mvpc_drw[[i_ind]]= test_adj(res_mvpc_drw, CPDAG)
  rp_ref[[i_ind]]= test_adj(res_mcar_pc,CPDAG)
  rp_pc[[i_ind]] = test_adj(res_pc,CPDAG)
  rp_td_pc[[i_ind]]=  test_adj(res_td_pc,CPDAG)
  
}

td_pc[[count]] = mean(shd_td_pc)
mvpc_drw[[count]] = mean(shd_mvpc_drw)
mvpc_permc[[count]] = mean(shd_mvpc_permc)
ref[[count]] = mean(shd_ref)
pc[[count]] = mean(shd_pc)

std_td_pc[[count]] = sd(shd_td_pc)
std_mvpc_drw[[count]] = sd(shd_mvpc_drw)
std_mvpc_permc[[count]] = sd(shd_mvpc_permc)
std_ref[[count]] = sd(shd_ref)
std_pc[[count]] = sd(shd_pc)

r_td_pc[[count]] = compute_rp(rp_td_pc)
r_mvpc_drw[[count]] = compute_rp(rp_mvpc_permc)
r_mvpc_permc[[count]] = compute_rp(rp_mvpc_drw)
r_ref[[count]] = compute_rp(rp_ref)
r_pc[[count]] = compute_rp(rp_pc)
count = count + 1

cat(paste(mean(shd_pc),sd(shd_pc),"\n"),file="output_mnar.txt",append=TRUE)
cat(paste(mean(shd_ref),sd(shd_ref),"\n"),file="output_mnar.txt",append=TRUE)
cat(paste(mean(shd_mvpc_permc),sd(shd_mvpc_permc),"\n"),file="output_mnar.txt",append=TRUE)
cat(paste(mean(shd_mvpc_drw),sd(shd_mvpc_drw),"\n"),file="output_mnar.txt",append=TRUE)
cat(paste(mean(shd_td_pc),sd(shd_td_pc),"\n"),file="output_mnar.txt",append=TRUE)

cat(paste(compute_f1(rp_pc),"\n"),file="output_mnar.txt",append=TRUE)
cat(paste(compute_f1(rp_ref),"\n"),file="output_mnar.txt",append=TRUE)
cat(paste(compute_f1(rp_mvpc_permc),"\n"),file="output_mnar.txt",append=TRUE)
cat(paste(compute_f1(rp_mvpc_drw),"\n"),file="output_mnar.txt",append=TRUE)
cat(paste(compute_f1(rp_td_pc),"\n"),file="output_mnar.txt",append=TRUE)
