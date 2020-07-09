library(pcalg)
library(mvtnorm)
library(weights)
library(ks)
library(e1071)  # For generating all combinations of a binary vector
library(weights)  # For kernel density estimate
library(ks)

library(DescTools)
library(mipfp)
#******************Functions for Synthetic Data Generation******************
gen_data <- function(num_sample,
                     mode='mar',
                     num_var=20, 
                     num_extra_e=3, 
                     num_m = 6, 
                     seed=1000,
                     p_missing_h=0.9, 
                     p_missing_l=0.1){
  # p: number of variables 
  # n: data sample size
  # mode: different methods to generate data sets with different missingness mechanisms, such as MCAR, MAR and MNAR
  # num_extra_e: number of the variables with missing values that lead to wrong results
  # num_m: number of variables with missing values
  # seed: random seeds
  # p_missing_h: The probability of missing values when the missingness condition is satisfied, e.g., missingness indicator R = 1
  # p_missing_l: The probability of missing values when the missingness condition is not satisfied, e.g., missingness indicator R = 0
  # **************
  # Return: 
  # data_complete: the complete data set
  # data_m: the data set containing missing values that generated with mode function
  # data_ref: the MCAR data set as reference 
  # ground_truth: the ground-true DAG, CPDAG, 
  
  set.seed(seed) # one seed has a corresponding random graph (seed controls the graph)
  myDAG <- randomDAG(num_var,2/(num_var-1))
  # Make sure the collider is the cause of R-variable of its parent variable
  data_del <-gen_del(num_sample, myDAG, mode, 
                     num_m=num_m, 
                     num_var=num_var, 
                     num_extra_e =num_extra_e,
                     p_missing_h, p_missing_l)
  
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

gen_m_prt <- function(DAG, mode='mar',
                      num_var=20, 
                      num_extra_e=3, 
                      num_m = 6){
  # Given a DAG, return a list of missingness indicators and their parents
  
  cldr <- detect_colliders(DAG)
  cldr_prt <- detect_colliders_prt(DAG, cldr)
  # Choose missingness inidcator and their parents
  if(mode=='mar'){
    p_m <- create_mar_ind(cldr,cldr_prt,num_var, num_extra_e, num_m)
  }
  else if(mode == 'mnar'){
    p_m <- create_mnar_ind(cldr,cldr_prt,num_var, num_extra_e, num_m)
  }else{
    p_m <- create_mul_mar_ind(cldr,cldr_prt,num_var, num_extra_e, num_m)
  }
  
  ms = p_m$ms
  prt_ms = p_m$prt_ms
  return(list(ms = ms, prt_ms=prt_ms))
}

gen_del <- function(n,myDAG,mode='mar',
                    num_m=6,
                    num_var=20, 
                    num_extra_e=3,
                    p_missing_h=0.9, p_missing_l=0.01){
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
  
  num_samples <- n
  prt_m = gen_m_prt(myDAG,mode,
                    num_m=num_m, 
                    num_var=num_var, 
                    num_extra_e =num_extra_e )
  
  
  # 3. Generate missing values according the the values of missingness indicators
  data_com_m_ref = generate_missing_values(num_var,n,myDAG, prt_m$ms, prt_m$prt_ms,
                                           p_missing_h, p_missing_l)
  data_m = data_com_m_ref$data_m
  data_ref = data_com_m_ref$data_ref
  data_complete = data_com_m_ref$data_complete
  
  return(list(data_complete=data_complete,
              data_m=data_m,
              data_ref=data_ref,
              parent_m_ind=prt_m$prt_ms,
              m_ind=prt_m$ms)) 
}

generate_missing_values <- function(p,n,myDAG,m_ind,parent_m_ind, p_missing_h=0.9, p_missing_l=0.01){
  # Give the parents of missingness indicators, and the missingness indcators
  # The missing values are generated whn the parent values are in the bottom XX percentage.
  data <- rmvnorm(n, mean=rep(0,p), sigma=trueCov(myDAG))
  data_m=data
  data_mcar=data
  for(i in c(1:length(m_ind))){
    # Choose lower "bottom_p" percentage of the values
    bottom_p <- runif(1, min = 0.1, max = 0.7)
    r <- mis_cal_ind(data[,parent_m_ind[i]], qnorm(bottom_p), p_missing_h, p_missing_l)
    data_m[r, m_ind[i]] = NA
    r_mcar<-sample(r)
    data_mcar[r_mcar, m_ind[i]] = NA
  }
  return(list(data_complete=data, 
              data_m=data_m,
              data_ref=data_mcar))
}

mis_cal_ind<-function(x, bottom_p, p_missing_h=0.9, p_missing_l=0.01){
  ind = x < bottom_p
  h_x = rbinom(length(x),1,p_missing_h) ==1
  l_x = rbinom(length(x),1,p_missing_l) ==1
  out = l_x
  out[ind] = h_x[ind]
  out
} 

detect_colliders <- function(myDAG){
  m <- as(myDAG,"matrix")
  m[m>0]<-1
  ind <- colSums(m,1)>1
  num_vars<-length(m[,1])
  x <- c(1:num_vars)
  x[ind]
}

detect_colliders_prt <- function(DAG, cldr){
  m <- as(DAG,"matrix")
  a = 1:length(m[,1])
  cldr_prt = list()
  count = 1
  for(i in cldr){
    cldr_prt[[count]]=a[m[,i]>0] # contain the situation where  "== 1"
    count = count + 1 
  }
  cldr_prt
}

prtm.groundtrue <- function(m,prt){
  # get ground-true missingess parents
  # return prt_m as a data.frame for mvpc
  # Example: ground-true missingess parents
  # prt = gen_result_list$ground_truth$parent_m_ind
  # m = gen_result_list$ground_truth$m_ind
  # prtm.groundtrue(m,prt)
  
  prt_m<-data.frame(m=m)
  prt1 = list()
  for(i in 1:length(prt)){
    prt1[[i]] = c(prt[i])
  }
  prt_m[['prt']]<-prt1
  return(prt_m)
}

create_mnar_ind <- function(cldr,cldr_prt,num_var=20, num_extra_e=3, num_m = 6){
  prt_ms = c()
  ms = c()
  count = 1
  ## Create MAR 
  for(i in 1:length(cldr)){
    for(pr in cldr_prt[[i]]){
      if(length(ms) == 0){
        ms[count] = pr
        prt_ms[count] = cldr[i]
        count = count + 1
      }
      else{
        if((!(pr %in% prt_ms)) && (! (pr %in% ms))){
          ms[count] = pr
          prt_ms[count] = cldr[i]
          count = count + 1
        }
      }
    }
  }
  
  # Only involve "num_extra_e" number of colliders
  if(count > (num_extra_e+1)){
    ind_cld = sample(1:length(ms))
    ms = ms[ind_cld[1:num_extra_e]]
    prt_ms = prt_ms[ind_cld[1:num_extra_e]]
    count = (num_extra_e+1)
  }
  
  # Add MCAR over MAR for generating MNAR
  
  ind_rd = count
  # Start from missingness indicators
  left_ind_prt = setdiff(1:num_var, prt_ms) # The parent can be a collider again.
  left_ind_prt = setdiff(left_ind_prt, ms)   # MNAR: one parent of a missingness indictor
  left_ind_prt = sample(left_ind_prt) 
  end_for = num_m - length(ms)
  
  countm = count
  for(i in 1:end_for){
    # Append prt not MAR, i.e., not collider
    if(!(i > num_extra_e)){ms[countm] = prt_ms[i]}
    else{ms[countm] = left_ind_prt[i]}
    countm = countm + 1
  }
  
  countp = count
  left_ind_prt = setdiff(1:num_var, cldr)
  for(i in 1:end_for){
    # Append missingness indicator not self-masking
    not_self_masking = sample(setdiff(left_ind_prt, ms[countp]))[1]
    prt_ms[countp] = not_self_masking
    countp = countp + 1
  }
  
  return(list(ms = ms, prt_ms=prt_ms))   
}

create_mar_ind <- function(cldr,cldr_prt,num_var=20, num_extra_e=3, num_m = 6){
  prt_ms = c()
  ms = c()
  count = 1
  
  for(i in 1:length(cldr)){
    for(pr in cldr_prt[[i]]){
      if(length(ms) == 0){
        ms[count] = pr
        prt_ms[count] = cldr[i]
        count = count + 1
      }
      else{
        if((!(pr %in% prt_ms)) && (! (pr %in% ms))){
          ms[count] = pr
          prt_ms[count] = cldr[i]
          count = count + 1
        }
      }
    }
  }
  
  # Only involve "num_extra_e" number of colliders
  if(count > (num_extra_e+1)){
    ind_cld = sample(1:length(ms))
    ms = ms[ind_cld[1:num_extra_e]]
    prt_ms = prt_ms[ind_cld[1:num_extra_e]]
    count = (num_extra_e+1)
  }
  
  ind_rd = count
  left_ind_prt = 1:num_var  # The parent can be a collider again.
  left_ind_prt = setdiff(left_ind_prt, ms)
  left_ind_prt = setdiff(left_ind_prt, prt_ms)
  left_ind_prt = sample(left_ind_prt) # used for missingness indicators -- not collider, not in the ms and prt_ms
  
  end_for = num_m - length(ms)
  countp = count
  for(i in 1:end_for){
    # Append prt not MAR, i.e., not collider
    prt_ms[countp] = left_ind_prt[i]
    countp = countp + 1
  }
  
  left_ind_m = setdiff(1:num_var, ms) 
  left_ind_m = setdiff(left_ind_m, prt_ms)
  left_ind_m = sample(left_ind_m)
  countm = count
  
  for(i in 1:end_for){
    # Append prt not MAR, i.e., not collider
    ms[countm] = left_ind_m[i]
    countm = countm + 1
  }
  
  return(list(ms = ms, prt_ms=prt_ms))   
  
}

create_mul_mar_ind <- function(cldr,cldr_prt,num_var=20, num_extra_e=3, num_m = 6){
  prt_ms = c()
  ms = c()
  count = 1
  
  for(i in 1:length(cldr)){
    for(pr in cldr_prt[[i]]){
      if(length(ms) == 0){
        ms[count] = pr
        prt_ms[count] = cldr[i]
        count = count + 1
      }
      else{
        if((!(pr %in% prt_ms)) && (! (pr %in% ms))){
          ms[count] = pr
          prt_ms[count] = cldr[i]
          count = count + 1
        }
      }
    }
  }
  
  # Only involve "num_extra_e" number of colliders
  if(count > (num_extra_e+1)){
    ind_cld = sample(1:length(ms))
    ms = ms[ind_cld[1:num_extra_e]]
    prt_ms = prt_ms[ind_cld[1:num_extra_e]]
    count = (num_extra_e+1)
  }
  
  ind_rd = count
  left_ind_prt = 1:num_var  # The parent can be a collider again.
  left_ind_prt = setdiff(left_ind_prt, ms)
  left_ind_prt = setdiff(left_ind_prt, prt_ms)
  left_ind_prt = sample(left_ind_prt) # used for missingness indicators -- not collider, not in the ms and prt_ms
  
  end_for = num_m - length(ms)
  countp = count
  for(i in 1:end_for){
    # Append prt not MAR, i.e., not collider
    prt_ms[countp] = left_ind_prt[i]
    countp = countp + 1
  }
  
  left_ind_m = setdiff(1:num_var, ms) 
  left_ind_m = setdiff(left_ind_m, prt_ms)
  left_ind_m = sample(left_ind_m)
  countm = count
  
  for(i in 1:end_for){
    # Append prt not MAR, i.e., not collider
    if(i < count){
      ms[countm] = ms[i]  # repeat the missingness indicator to generate multiple cause
    }
    else{
      ms[countm] = left_ind_m[i]    
    }
    countm = countm + 1
  }
  
  return(list(ms = ms, prt_ms=prt_ms))   
  #  [m1, m2, m3, m1, m2, m3] 
  #  [prt[1], prt[2], prt[3],prt[4], prt[5], prt[6]]
}

load_bin_data<-function(fdata){
  read.table(fdata, header=TRUE, sep="\t", stringsAsFactors = FALSE) 
}

load_bin_graph<-function(graph.file){
  res <- readLines(graph.file)
  num.var <- length(strsplit(res[2], ",")[[1]])
  graph <- matrix(0L, nrow = num.var, ncol = num.var)
  for(row in res){
    if(grepl("->", row, fixed=TRUE)){
      itms <- strsplit(row, " ")[[1]]
      cause = as.integer(gsub("[^0-9.]", "",  itms[2]))
      effect = as.integer(gsub("[^0-9.]", "",  itms[4]))
      graph[cause,effect] = 1
    }
  }
  DAG <- as(graph,'graphNEL')
}