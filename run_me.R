# The path of "R-proj"
proj_path<-'/Users/ruibo/Desktop/mvpc/mvpc_v2/MVPC'
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")

source(paste(src_path,'/simulator.R',sep=""))
gen_result_list<-gen_data(20,10,"mnar")

data_m = gen_result_list$data_m

# Detecting variables with missing values in the data set
R<-get_m_ind(data_m)

# Detect the parents of the variables with missing values
sup_var<-get_prt_m_ind(R, data=data_m)
