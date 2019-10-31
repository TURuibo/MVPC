proj_path<-getwd()
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")

source(paste(src_path,'/all_functions.R',sep=""))

set.seed(777)
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
nrep = 5
for(i_g in 1:10){
  print(paste('graph=', i_g))
  data.name = paste('data',i_g,sep='')
  graph.name = paste('graph',i_g,sep='')
  
  data.file <- paste(data_path,'/',data.name,'.txt',sep="")
  graph.file<- paste(data_path,'/',graph.name,'.txt',sep="")
  
  data = load_bin_data(data.file)
  data = data[1:10000,]
  DAG = load_bin_graph(graph.file)
  CPDAG = dag2cpdag(DAG)
  
  for(i_exp in 1:nrep){
    print(paste('rep=', i_exp))
    # Detect colliders and  Colliders' parents
    cldr <- detect_colliders(DAG)
    cldr_prt <- detect_colliders_prt(DAG, cldr)
    # Choose missingness inidcator and their parents

    # p_m <- create_mar_ind(cldr,cldr_prt)
    p_m <- create_mnar_ind(cldr,cldr_prt)
    ms = p_m$ms
    prt_ms = p_m$prt_ms
    
    # Generate missing values 
    mask = data != data
    
    for(i in 1:length(ms)){
      nsample = nrow(data)
      p_1 = runif(1,0.1,0.5)
      m_ind = rbinom(nsample, 1, p_1)==1
      m_ind = (data[,prt_ms[i]] == 1) & m_ind
      mask[,ms[i]] = m_ind
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

    i_ind = (i_g-1)*nrep + i_exp
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
}

mean(shd_td_pc)
mean(shd_mvpc_drw)
mean(shd_mvpc_permc)
mean(shd_ref)
mean(shd_pc)

sd(shd_td_pc)
sd(shd_mvpc_drw)
sd(shd_mvpc_permc)
sd(shd_ref)
sd(shd_pc)

compute_rp(rp_td_pc)
compute_rp(rp_mvpc_permc)
compute_rp(rp_mvpc_drw)
compute_rp(rp_ref)
compute_rp(rp_pc)

