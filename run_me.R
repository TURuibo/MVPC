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

dag_lw <- pc(suffStat_ref_lw_del, gaussCItest_tw_del, alpha=0.01, p=num_var)
shd_ref_lw = shd(dag_lw, myCPDAG)

dag_tw <- pc(suffStat_ref, gaussCItest_tw_del, alpha=0.01, p=num_var)
shd_ref_tw = shd(dag_tw, myCPDAG)

suffStat_m <- list(data=data_m)
suffStat_m_lw_del<- list(data=test_wise_deletion(1:num_var, data_m))

dag_m_lw <- pc(suffStat_m_lw_del, gaussCItest_tw_del, alpha=0.01, p=num_var)
dag_m_tw <- pc(suffStat_m, gaussCItest_tw_del, alpha=0.01, p=num_var)
shd_m_lw = shd(dag_m_lw, myCPDAG)
shd_m_tw = shd(dag_m_tw, myCPDAG)

print(shd_comp)
print(shd_ref_lw)
print(shd_ref_tw)
print(shd_m_lw)
print(shd_m_tw)

# ********* Missing Value PC (MVPC) *********
mvpc<-function(suffStat, indepTest, alpha, labels, p, 
               fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE, m.max = Inf, 
               u2pd = c("relaxed", "rand", "retry"), 
               skel.method = c("stable", "original", "stable.fast"), 
               conservative = FALSE, maj.rule = FALSE, solve.confl = FALSE, numCores = 1, verbose = FALSE){
  # MVPC step1: Test-wise skeleton search.
  # (to extract the necessary information for correcting the wrong edges in the result and initialize the PC algorithm).
  
  
  # MVPC step2: Detect parents of missingness indicators.
  sup_var<-get_prt_m_ind(data=suffStat$data) # "data" is "data_m" which containing missing values.
  
  # MVPC step3: 
  # a) Run PC algorithm with the 1st step skeleton; 
  # b) Correct the wrong edges of it with permutation-based CI test (PermCCItest) and density ratio weighted CI test (DRWCItest).
  
}




