# The path of "R-proj"
proj_path<-'/Users/ruibo/Desktop/mvpc/mvpc_v2/MVPC'
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")

source(paste(src_path,'/all_functions.R',sep=""))

# ********* Synthethic data generation ********* 
num_var = 20
num_sample = 100000
rdm_seed = 10

gen_result_list<-gen_data(num_var,num_sample,"mnar",rdm_seed)

data_complete = gen_result_list$data_complete
data_m = gen_result_list$data_m
data_ref = gen_result_list$data_ref
myCPDAG = gen_result_list$ground_truth$cpdag
collider = gen_result_list$ground_truth$collider
prt = gen_result_list$ground_truth$parent_m_ind

print(collider)
print(prt) 
print(gen_result_list$ground_truth$m_ind)

# ********* Test-wise deletion PC ********* 
suffStat_complete <- list(data=data_complete)
dag <- pc(suffStat_complete, gaussCItest_tw_del, alpha=0.01, p=num_var)
shd_comp = shd(dag, myCPDAG)

suffStat_ref <- list(data=data_ref)
suffStat_ref_lw_del<- list(data=test_wise_deletion(1:num_var, data_ref))
dag <- pc(suffStat_ref_lw_del, gaussCItest_tw_del, alpha=0.01, p=num_var)
shd_ref = shd(dag, myCPDAG)

suffStat_m <- list(data=data_m)
suffStat_m_lw_del<- list(data=test_wise_deletion(1:num_var, data_m))
dag <- pc(suffStat_m_lw_del, gaussCItest_tw_del, alpha=0.01, p=num_var)
shd_m = shd(dag, myCPDAG)

print(shd_comp)
print(shd_ref)
print(shd_m)


# ********* Missing Value PC (MVPC) ********* 
# Detecting variables with missing values in the data set
R<-get_m_ind(data_m)

# Detect the parents of the variables with missing values
sup_var<-get_prt_m_ind(R, data=data_m)
