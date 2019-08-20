# The path of "R-proj"
proj_path<-getwd()
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")

source(paste(src_path,'/all_functions.R',sep=""))

# ********* Synthethic Binary Data Generation ********* 
data = load_bin_data()
graph = load_bin_graph()

# ********* Missingness generation ********* 

# ********* Test-Wise Deletion ********* 
res_pc<-pc(suffStat_m, binCItest_td, alpha=0.01, p=num_var)

# ********* Correction ********* 

res_mvpc<-mvpc(suffStat_m, binCItest_td, binPermCCItest, alpha=0.01, p=num_var)

