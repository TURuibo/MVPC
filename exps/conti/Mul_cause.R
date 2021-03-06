proj_path<-getwd()
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")

source(paste(src_path,'/CITest.R',sep=""))
source(paste(src_path,'/Evaluation.R',sep=""))
source(paste(src_path,'/MissingValuePC.R',sep=""))
source(paste(src_path,'/SyntheticDataGeneration.R',sep=""))

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

## ********* Synthethic Binary Data Generation ********* 
num_var=20
num_sample = 50000
count = 1

print(paste('num_sample: ', num_sample))
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
for(graph_ind in 1:10){
  print(graph_ind)
  gen_result_list<-gen_data(num_sample = num_sample, 
                            mode = "mul_mar",
                            num_var=num_var, 
                            num_extra_e=5, 
                            num_m = 10, 
                            seed = graph_ind)
  
  data_complete = gen_result_list$data_complete
  data_m = gen_result_list$data_m
  myCPDAG = gen_result_list$ground_truth$cpdag
  prt = gen_result_list$ground_truth$parent_m_ind
  m = gen_result_list$ground_truth$m_ind
  m_u = unique(m)
  prt_m <- data.frame(m=m_u)
  prt_m['prt'] = list()
  
  for(i in 1:length(m_u)){
    prt_m[['prt']][[i]] = c(prt[m==m_u[i]])
  }
  
  suffStat_m <- list(data=data_m, prt_m=prt_m)
  suffStat  = list(C = cor(data_complete),n=num_sample)
  
  ## ********* PermC Correction *********
  res.mvpc.permc <-mvpc(suffStat_m,
                        gaussCItest.td, 
                        gaussCItest.permc,
                        prt_m=prt_m, alpha=0.01, p=num_var)
  
  ## ********* DRW Correction *********
  res.mvpc.drw <-mvpc(suffStat_m,
                      gaussCItest.td,
                      gaussCItest.drw,
                      prt_m=prt_m, alpha=0.01, p=num_var)
  
  ## ********* Complete data evaluation *********
  res.comp<-pc(suffStat, gaussCItest, alpha=0.01, p=num_var)
  
  ## ********* Test-Wise Deletion *********
  sample_size <<- c()
  suffStat_tw = list(data=data_m)
  res.td<-pc(suffStat_tw, gaussCItest.td.ref, alpha=0.01, p=num_var)
  sample_size <- floor(mean(sample_size))
  
  ## ********* MCAR Complete data evaluation *********
  suffStat_ref  = list(C = cor(data_complete[1:sample_size,]),n=sample_size)
  res.ref<-pc(suffStat_ref, gaussCItest, alpha=0.01, p=num_var)
  
  i_ind = graph_ind
  
  shd_mvpc_permc[i_ind] =shd(res.mvpc.permc,myCPDAG)
  shd_mvpc_drw[i_ind] =shd(res.mvpc.drw,myCPDAG)
  shd_ref[i_ind] =shd(res.ref,myCPDAG)
  shd_pc[i_ind] = shd(res.comp,myCPDAG)
  shd_td_pc[i_ind] =shd(res.td,myCPDAG)
  
  rp_mvpc_permc[[i_ind]]= test_adj(res.mvpc.permc, myCPDAG)
  rp_mvpc_drw[[i_ind]]= test_adj(res.mvpc.drw, myCPDAG)
  rp_ref[[i_ind]]= test_adj(res.ref,myCPDAG)
  rp_pc[[i_ind]] = test_adj(res.comp,myCPDAG)
  rp_td_pc[[i_ind]]=  test_adj(res.td,myCPDAG)
}
print(c(mean(shd_pc),sd(shd_pc)))
print(c(mean(shd_ref),sd(shd_ref)))
print(c(mean(shd_mvpc_permc),sd(shd_mvpc_permc)))
print(c(mean(shd_mvpc_drw),sd(shd_mvpc_permc)))
print(c(mean(shd_td_pc),sd(shd_td_pc)))

print(compute_f1(rp_pc))
print(compute_f1(rp_ref))
print(compute_f1(rp_mvpc_permc))
print(compute_f1(rp_mvpc_drw))
print(compute_f1(rp_td_pc))

