# The path of "R-proj"
proj_path<-getwd()
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")

source(paste(src_path,'/all_functions.R',sep=""))

# ********* Synthethic data generation ********* 

num_var=20
num_sample = 100000
gen_result_list<-gen_data(num_sample = num_sample, 
                          mode = "mnar",
                          num_var=num_var, 
                          num_extra_e=5, 
                          num_m = 10, 
                          seed = 777)

data_complete = gen_result_list$data_complete
data_m = gen_result_list$data_m
myCPDAG = gen_result_list$ground_truth$cpdag

# ground-true missingess parents
prt = gen_result_list$ground_truth$parent_m_ind
m = gen_result_list$ground_truth$m_ind
prtm.groundtrue(m,prt)

# ********* MVPC *********

suffStat_m <- list(data=data_m)
suffStat  = list(C = cor(data_complete),n=num_sample)

res_tw<-pc(suffStat_m,
           gaussCItest.td,
           alpha=0.01, 
           p=num_var)

res_com_pc<-pc(suffStat, 
               gaussCItest, 
               alpha=0.01, 
               p=num_var)

res.mvpc.permc <-mvpc(suffStat_m,
                      gaussCItest.td, 
                      gaussCItest.PermC,
                      alpha=0.01,
                      p=num_var)

res.mvpc.drw <-mvpc(suffStat_m,
                    gaussCItest.td,
                    gaussCItest.drw,
                    alpha=0.01, 
                    p=num_var)


shd(res.mvpc.drw,myCPDAG)
shd(res.mvpc.permc,myCPDAG)
shd(res_tw,myCPDAG)
shd(res_com_pc,myCPDAG)
