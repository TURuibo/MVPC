library(pcalg)
library(Rfast)
load_sim_data <- function(p,n,mode=gen_mar){
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
  set.seed(1000)
  data_complete <- rmvnorm(n, mu=rep(0,p), sigma=trueCov(myDAG))
  
  # Make sure the collider is the cause of R-variable of its parent variable
  data_del <-mode(data_complete,myDAG)
  data_m = data_del$data
  # Also provide reference, a MCAR data set and the ground-truth (DAG, missingness indicator, and parents of the missingness indicators ).
  data_ref <- data.frame(data_del$data_mcar)
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

gen_mar <- function(data, myDAG){
  return(list(data=NULL,
              data_ref=NULL,
              parent_m_ind=NULL,
              m_ind=NULL
              )) 
}

load_sim_data(20,10)

# gen_mar <- function(data, myDAG){
#   num_m<-sample(6:10,1)
#   size <- dim(data)
#   num_var <- size[2]
#   num_samples <- size[1]
#   m_ind <- sample(c(1:num_var))
#   colliders <- find_colliders(myDAG)
#   if(length(colliders)>num_m){
#     colliders <- colliders[1:num_m]
#   }
#   not_colliders <- setdiff(c(1:num_var),colliders)
#   
#   sup_ind<-rep(0,1)
#   mar_ind<-rep(0,1)
#   count = 1 
#   i=1
#   sep<-ceiling(length(colliders)/2)
#   
#   for(cldr in colliders){
#     if(i <= sep){prts <- find_parents(cldr, myDAG)}
#     else{
#       prts <- c(1:num_var)
#     }
#     i=i+1
#     prts<-setdiff(prts,mar_ind)
#     prts<-setdiff(prts,colliders)
#     if(length(prts)>0){
#       mar_ind[count] <- prts[1]
#       sup_ind[count] <- cldr
#       count<-count+1
#     }
#     else{
#       next
#     }
#   }
#   # what left
#   if(length(sup_ind)>=num_m){
#     sup_ind<-sup_ind[1:num_m]
#     mar_ind<-mar_ind[1:num_m]
#   }
#   else{
#     x = c(1:num_var)
#     x<-setdiff(x,mar_ind)
#     x<-setdiff(x,sup_ind)
#     if(length(x)>0){
#       x<-sample(x)
#       sup_ind[count:num_m] <- x[1:(num_m-count+1)]
#       mar_ind[count:num_m] <- x[(num_m-count+1+1):(2*(num_m-count+1))]
#     }
#   }
#   
#   data_mcar<-data
#   for(i in c(1:num_m)){
#     r_bottom <- runif(1, min = 0.1, max = 0.5)
#     r <- data[,sup_ind[i]] < qnorm(r_bottom)
#     
#     data[r==1, mar_ind[i]] = NA
#     r_mcar<-sample(r)
#     data_mcar[r_mcar==1, mar_ind[i]] = NA
#   }
#   return(list(data=data,
#               data_mcar=data_mcar,
#               parent_m_ind=sup_ind,
#               m_ind=mar_ind))
# }
