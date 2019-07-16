library(pcalg)
library(mvtnorm)

#******************Functions for synthetic data generation******************
gen_data <- function(p,n,mode='mar',seed=1000){
  # p: number of variables 
  # n: data sample size
  # mode: different methods to generate data sets with different missingness mechanisms, such as MCAR, MAR and MNAR
  # **************
  # Return: 
  # data_complete: the complete data set
  # data_m: the data set containing missing values that generated with mode function
  # data_ref: the MCAR data set as reference 
  # ground_truth: the ground-true DAG, CPDAG, 
  
  set.seed(seed) # one seed has a corresponding random graph (seed controls the graph)
  myDAG <- randomDAG(p,2/(p-1))
  ## Make sure the graph contains collider
  # while(!check_collider(myDAG)){
  #   myDAG <- randomDAG(p,2/(p-1))
  # }
  
  # Make sure the collider is the cause of R-variable of its parent variable
  data_del <-gen_del(p,n,myDAG,mode)
  data_m = data_del$data_m
  data_complete = data_del$data_complete
  data_ref <- data.frame(data_del$data_ref)
  # Also provide reference, a MCAR data set and the ground-truth (DAG, missingness indicator, and parents of the missingness indicators ).
  ground_truth <- list(
    dag=myDAG,
    cpdag=dag2cpdag(myDAG),
    collider=detect_colliders(myDAG),
    parent_m_ind=data_del$parent_m_ind,
    m_ind=data_del$m_ind
  )

  return(list(data_complete=data_complete,
              data_m=data_m,
              data_ref=data_ref,
              ground_truth=ground_truth))
}

gen_del <- function(p,n,myDAG,mode='mar',num_m=6){
  # Constraint of MAR:
  # 1. the missingness indicator is caused by some substantive variables;
  # 2. the parents of missingness indicator have no missing values in the data set.
  # Constraint of our implememntation:
  # Here we only choose a single parent for a missingness indicator.
  # ********************
  # Constraints of MNAR:
  # 1. no self-masking: X -> Rx
  # 2. parents of missingness indicators must also have missing values in the data set.
  # 3. the limitation of the correction method 
  # ********************
  
  num_var <- p
  num_samples <- n
  # 1. Choose the missingness indicators/ choose which variables contain missing values 
  m_ind = select_m_ind(num_var,num_m)
  
  # 2. Choose the cause of missingess indicators
  if(mode=="mar"){
    prt_m = select_prt_mar_ind(m_ind,num_var)
  }else{
    prt_m = select_prt_mnar_ind(m_ind,num_var)
  }
  
  # 3. Generate missing values according the the values of missingness indicators
  data_com_m_ref = generate_missing_values(p,n,myDAG,m_ind,prt_m)
  data_m = data_com_m_ref$data_m
  data_ref = data_com_m_ref$data_ref
  data_complete = data_com_m_ref$data_complete
  
  return(list(data_complete=data_complete,
              data_m=data_m,
              data_ref=data_ref,
              parent_m_ind=prt_m,
              m_ind=m_ind)) 
}

select_m_ind<-function(num_var,num_m){
  return(sample(x=1:num_var,size=num_m,replace=FALSE))
}

select_prt_mar_ind<-function(m_ind,num_var){
  return(sample(x=setdiff(1:num_var,m_ind),size=length(m_ind),replace=TRUE))
}

select_prt_mnar_ind<-function(m_ind,num_var){
  # No self-masking
  # only the variable contains missing values can be the parents of missingness indicators
  prts = c()
  for(i in m_ind){
    newelem<-sample(setdiff(m_ind,i),size=1)
    prts <- c(prts, newelem) 
  }
  return(prts)
}

generate_missing_values<-function(p,n,myDAG,m_ind,parent_m_ind){
  # Give the parents of missingness indicators, and the missingness indcators
  # The missing values are generated whn the parent values are in the bottom XX percentage.
  data <- rmvnorm(n, mean=rep(0,p), sigma=trueCov(myDAG))
  data_m=data
  data_mcar=data
  for(i in c(1:length(m_ind))){
    # Choose lower "bottom_p" percentage of the values
    bottom_p <- runif(1, min = 0.1, max = 0.5)
    r <- data[,parent_m_ind[i]] < qnorm(bottom_p)
    data_m[r==1, m_ind[i]] = NA
    r_mcar<-sample(r)
    data_mcar[r_mcar==1, m_ind[i]] = NA
  }
  return(list(data_complete=data, 
              data_m=data_m,
              data_ref=data_mcar))
}

detect_colliders<-function(myDAG){
  m <- as(myDAG,"matrix")
  m[m>0]<-1
  ind <- colSums(m,1)>1
  num_vars<-length(m[,1])
  x <- c(1:num_vars)
  x[ind]
}

# Example of generating synthetic data
# gen_data(20,10,"mnar")

#****************** Independence test ****************** 
gaussCItest_tw_del <- function(x, y, S, suffStat) {
  ## Conditional independence test between continuous variables with deletion methods
  ## test P(x,y|S)
  ## suffStat: the class contains the dataset and other necessary variables
  ##    suffStat$data: the dataset
  ##--------------
  ## Return: the p-value of the test 
  
  data = test_wise_deletion(c(x,y,S),suffStat$data)
  z <- zStat(x,y,S, 
             C = cor(data), n = length(data[,1]))
  2*pnorm(abs(z), lower.tail = FALSE)
}

test_wise_deletion <-function(var_ind, data){
  ## Delete the rows of given variables (var_ind) if there is a missing value in a row
  ## var_ind: variables in the current conditonal independence test
  ## data: the whole data set 
  ##--------------------------------
  ## Return: the deleted dataset of the variables in the CI test 
  
  not_del_ind = c(rep(TRUE,length(data[,1])))
  for (var in var_ind){
    if(anyNA(data[,var])){
      not_del_ind = not_del_ind & !is.na(data[,var])
    }
  }  
  return(data[not_del_ind,])
}

PermCCItest <- function(x, y, S, suffStat){}

DRWCItest <- function(x, y, S, suffStat){}

#****************** Functions of extracting necessary information for MVPC ******************
get_m_ind <- function(data){
  ## Check whether the current testing variables containing missing variables 
  ## Value: return TRUE or FALSE
  size <- dim(data)
  num_var <- size[2]
  m_ind = c()
  for (var in 1:num_var){
    if(anyNA(data[,var])){
      m_ind = c(m_ind,var)
    }
  }
  return(m_ind)
}

get_prt_m_ind<-function(data, alpha=0.01, mode=Anovatest){
  ## Detect causes of the missingness indicators
  ## R: 
  ## data:
  ## mode: Anova test to test the conditional independence
  ## ---------------------------
  ## Return: the index of parents of missingness indicators
  R<-get_m_ind(data_m)
  m=c()
  prt=list()
  suffStat = list(data=data)
  
  count=1
  for(R_ind in R){
    sup_test<- skeleton3(suffStat, indepTest=mode, R_ind=R_ind, alpha = alpha, verbose=FALSE)
    if(length(sup_test)!=0){
      m[count] <- R_ind
      prt[[count]]<-sup_test
      count = count + 1 
    }
  }
  prt_m<-data.frame(m=m)
  prt_m[['sup']]<-sup
  return(prt_m)
}
