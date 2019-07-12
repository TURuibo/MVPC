library(pcalg)
library(Rfast)

load_sim_data <- function(p,n,mode=gen_mar,seed=1000){
  # p: number of variables 
  # n: data sample size
  # mode: different methods to generate data sets with different missingness mechanisms, such as MCAR, MAR and MNAR
  # **************
  # return: 
  # data_complete: the complete data set
  # data_m: the data set containing missing values that generated with mode function
  # data_ref: the MCAR data set as reference 
  # ground_truth: the ground-true DAG, CPDAG, 
  
  myDAG <- randomDAG(p,2/(p-1))
  ## Make sure the graph contains collider
  # while(!check_collider(myDAG)){
  #   myDAG <- randomDAG(p,2/(p-1))
  # }
  set.seed(seed)
  data_complete <- rmvnorm(n, mu=rep(0,p), sigma=trueCov(myDAG))
  
  # Make sure the collider is the cause of R-variable of its parent variable
  data_del <-mode(data_complete,myDAG)
  data_m = data_del$data_mar
  # Also provide reference, a MCAR data set and the ground-truth (DAG, missingness indicator, and parents of the missingness indicators ).
  data_ref <- data.frame(data_del$data_ref)
  ground_truth <- list(
    dag=myDAG,
    cpdag=dag2cpdag(myDAG),
    parent_m_ind=data_del$parent_m_ind,
    m_ind=data_del$m_ind
  )

  return(list(data_complete=data_complete,
              data_m=data_m,
              data_ref=data_ref,
              ground_truth=ground_truth))
}

gen_mar <- function(data, myDAG, num_m=6){
  # Constraint of MAR:
  # 1. the missingness indicator is caused by some substantive variables;
  # 2. the parents of missingness indicator have no missing values in the data set.
  # Constraint of our implememntation:
  # Here we only choose a single parent for a missingness indicator.
  
  size <- dim(data)
  num_var <- size[2]
  num_samples <- size[1]
  # 1. Choose the missingness indicators/ choose which variables contain missing values 
  m_ind = select_missingness_var(num_var,num_m)
  
  # 2. Choose the cause of missingess indicators
  prt_m = select_parent_var(m_ind,num_var)
  
  # 3. Generate missing values according the the values of missingness indicators
  data_mar_ref = generate_missing_values(data,m_ind,prt_m)
  data_mar = data_mar_ref$data_mar
  data_ref = data_mar_ref$data_ref
  
  return(list(data_mar=data_mar,
              data_ref=data_ref,
              parent_m_ind=prt_m,
              m_ind=m_ind
              )) 
}

select_missingness_var<-function(num_var,num_m){
  return(sample(x=1:num_var,size=num_m,replace=FALSE))
}

select_parent_var<-function(m_ind,num_var){
  return(sample(x=setdiff(1:num_var,m_ind),size=length(m_ind),replace=TRUE))
}

generate_missing_values<-function(data,m_ind,parent_m_ind){
  data_mar=data
  data_mcar=data
  for(i in c(1:length(m_ind))){
    # Choose lower "bottom_p" percentage of the values
    bottom_p <- runif(1, min = 0.1, max = 0.5)
    r <- data[,parent_m_ind[i]] < qnorm(bottom_p)
    data_mar[r==1, m_ind[i]] = NA
    r_mcar<-sample(r)
    data_mcar[r_mcar==1, m_ind[i]] = NA
  }
  return(list(data_mar=data_mar,data_ref=data_mcar))
}

load_sim_data(20,10)
