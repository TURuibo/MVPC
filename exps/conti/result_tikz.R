proj_path<-getwd()
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")

source(paste(src_path,'/all_functions.R',sep=""))
sample.size= c(500,1000,5000,10000,50000,100000)
tab = read.table(paste(proj_path,'/exps/conti/asy_conti_mnar_mar_mcar',sep=""))
for(exp in 1:10 ){
  for(group.id in 1:6){
    if(exp<6){
      cat(paste("(",sample.size[group.id],",", tab[10*(group.id-1)+exp,1],")+-=(0,",tab[10*(group.id-1)+exp,2],")\n",sep=""),file="output.txt",append=TRUE)
    }
    else{
      cat(paste("(",sample.size[group.id],",", tab[10*(group.id-1)+exp,1],")+-=(0,",tab[10*(group.id-1)+exp,2],")\n",sep=""),file="output.txt",append=TRUE)
    }
  }
}