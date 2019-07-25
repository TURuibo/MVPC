# The path of "R-proj"
proj_path<-'/Users/ruibo/Desktop/mvpc/mvpc_v2/MVPC'
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")

source(paste(src_path,'/all_functions.R',sep=""))

# ********* Synthethic data generation ********* 

num_var = 20
num_sample = 1000
rdm_seed = 1

gen_result_list<-gen_data(num_var,num_sample,"mnar",rdm_seed)

data_complete = gen_result_list$data_complete
data_m = gen_result_list$data_m
data_ref = gen_result_list$data_ref
myCPDAG = gen_result_list$ground_truth$cpdag
collider = gen_result_list$ground_truth$collider
prt = gen_result_list$ground_truth$parent_m_ind
m = gen_result_list$ground_truth$m_ind
gth = data.frame(m=m,prt=prt)
gth <- gth[order(gth$m),] 
cat("Number of colliders that are parents of missingness inidcators: ", length(intersect(prt, collider)))


# ********* MVPC *********

suffStat_m <- list(data=data_m)
res<-mvpc(suffStat_m, gaussCItest_td, alpha=0.01, p=num_var)
print(gth)
print(res)
res <- data.frame(res)
cat("The wrong number of causes: ", eva.detection(gth$prt,res$prt))
# ********* MVPC *********
