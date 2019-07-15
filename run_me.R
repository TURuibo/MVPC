# The path of "R-proj"
proj_path<-'/Users/ruibo/Desktop/mvpc/mvpc_v2/MVPC'
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")

source(paste(src_path,'/all_functions.R',sep=""))

# ********* Synthethic data generation ********* 
num_var = 10
num_sample = 10000
rdm_seed = 100

gen_result_list<-gen_data(num_var,num_sample,"mar",rdm_seed)

data_complete = gen_result_list$data_complete
data_m = gen_result_list$data_m
data_ref = gen_result_list$data_ref
myCPDAG = gen_result_list$ground_truth$cpdag

# ********* Test-wise deletion PC ********* 
suffStat <- list(data=data_complete)
dag <- pc(suffStat, gaussCItest_tw_del, alpha=0.01, p=num_var)
shd_comp = shd(dag, myCPDAG)

suffStat <- list(data=data_ref)
dag <- pc(suffStat, gaussCItest_tw_del, alpha=0.01, p=num_var)
shd_ref = shd(dag, myCPDAG)

suffStat <- list(data=data_m)
dag <- pc(suffStat, gaussCItest_tw_del, alpha=0.01, p=num_var)
shd_m = shd(dag, myCPDAG)

print(shd_comp)
print(shd_ref)
print(shd_m)

# ********* Missing Value PC (MVPC) ********* 
# Detecting variables with missing values in the data set
R<-get_m_ind(data_m)

# Detect the parents of the variables with missing values
sup_var<-get_prt_m_ind(R, data=data_m)
