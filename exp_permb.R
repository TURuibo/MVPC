# The path of "R-proj"
proj_path<-'/Users/ruibo/Desktop/mvpc/mvpc_v2/MVPC'
src_path<-paste(proj_path,'/src',sep="")
res_path<-paste(proj_path,'/result',sep="")
data_path<-paste(proj_path,'/data',sep="")

source(paste(src_path,'/all_functions.R',sep=""))

n = 100000
a =  rnorm(n, 0, 1)
b =  rnorm(n, 0, 1)
e =  rnorm(n, 0, 1)

data_t = data.frame(a)
data_t[,2]=b

cor(data_t)
suffStat_t  = list(C=cor(data_t), n=n)
gaussCItest(1, 2, c(), suffStat_t)

# Generating missing data
w  = a + b + e 
r = w >  0
data_t[,3]=w

data_t2 = data_t
data_t2[!r,1] = NA

suffStat_t2 = list(data=data_t2)
gaussCItest_td(1, 2, c(), suffStat_t2)

## Permutation-based correction methods 
# Step 1: Learning generaive model for {X, Y, S} to impute X, Y, and S
data <- test_wise_deletion(c(1, 2, 3), data_t2)
xnam <- paste0("data[,", 3,"]")
fmla <- as.formula(paste("data[,1] ~ 0 + ", paste(xnam, collapse= "+")))
fit <- lm(formula = fmla, data = data)
res <- residuals(fit) # residuals

fmla2 <- as.formula(paste("data[,2] ~ 0 + ", paste(xnam, collapse= "+")))
fit2 <- lm(formula = fmla2, data = data)
res2 <- residuals(fit2) # residuals

# Step 2: Shuffle the source "W" -- parents of the missingness indicators
w_p = sample(w)
w_p = w_p[1:length(res)]
data[,3] = w_p
# Step 3: Generate the virtual data follows the full data distribution P(X, Y, S)

vir = predict(fit,list(data)) + res
vir2 = predict(fit2,list(data)) + res2

data_t3 = data.frame(vir)
data_t3[,2]=vir2

suffStat_t3 = list(data=data_t3)
gaussCItest_td(1, 2, c(), suffStat_t3)

lr <- list()
lr[[1]] <- fit
# par(mfrow=c(1,2))
# p_data = data_t2
# plot(p_data[,1],p_data[,2], pch=19)
# p_data = data_t3
# plot(p_data[,1],p_data[,2], pch=19)
# 

gen_vir <- function(t_var,index,suffStat){
  ## with the linear regression model, predict the virtual data
  ## t_var: variables in the current conditional independence test
  ## index: concatenate(t_var, causes (not in t_var) of missingness indicators of t_var)
  ## ----------------------------------------
  ## Return: the virtual data of variables in the current conditional independence test
  
  ind = setdiff(c(1:length(index)),c(2))
  
  ## Original deleted indicator
  ori_del_ind <- find_del_ind(1:length(index),suffStat$data[index])
  
  data <- deletion(index,suffStat)
  xnam <- paste0("data[,", 3:length(index),"]")
  fmla <- as.formula(paste("data[,1] ~ 0 + ", paste(xnam, collapse= "+")))
  fit <- lm(formula = fmla, data = data[,ind])
  res <- residuals(fit) # residuals
  
  ind = setdiff(c(1:length(index)),c(1))
  xnam2 <- paste0("data[,", 3:length(index),"]")
  fmla2 <- as.formula(paste("data[,2] ~ 0 + ", paste(xnam, collapse= "+")))
  fit2 <- lm(formula = fmla2, data = data[,ind])
  res2 <- residuals(fit2) # residuals
  
  # shuffle all the entries of support variables 
  # and delete according to original R-var
  data_ = suffStat$data[,index]
  data = shuffle(data_,(length(t_var)+1):length(index))
  
  ## Virtual deleted indicator
  vir_del_ind <- find_del_ind2((length(t_var)+1):length(index),data)
  
  # data <- data[left_var,]  # delete later and data with orginal number of samples
  # Make sure the NA can be put at the predictor?
  num_sample <- length(data[,1])
  res_<-rep(0,num_sample)
  res_[ori_del_ind]<-res
  res2_<-rep(0,num_sample)
  res2_[ori_del_ind]<-res2
  
  vir = predict(fit,list(data[,3:length(index)])) + res_
  vir2 = predict(fit2,list(data[,3:length(index)])) + res2_
  
  data_ = suffStat$data[,t_var]
  data_[,1]=vir
  data_[,2]=vir2
  data_ = data_[vir_del_ind&ori_del_ind,]
  return(data_)
}
