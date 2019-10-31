library(R.matlab)
library(pcalg)
library(mvtnorm)
library(weights)
library(ks)
library(e1071)  # For generating all combinations of a binary vector

library(DescTools)
library(mipfp)
#******************Functions for Synthetic Data Generation******************
gen_m_prt <- function(DAG, mode='mar',
                      num_var=20, num_extra_e=3, num_m = 6){
  # Given a DAG, return a list of missingness indicators and their parents
  
  cldr <- detect_colliders(DAG)
  cldr_prt <- detect_colliders_prt(DAG, cldr)
  # Choose missingness inidcator and their parents
  if(mode=='mar'){
    p_m <- create_mar_ind(cldr,cldr_prt,num_var, num_extra_e, num_m)
  }
  else{
    p_m <- create_mnar_ind(cldr,cldr_prt,num_var, num_extra_e, num_m)
  }

  ms = p_m$ms
  prt_ms = p_m$prt_ms
  return(list(ms = ms, prt_ms=prt_ms))
}

gen_data <- function(p,n,mode='mar',
                     num_var=20, num_extra_e=3, num_m = 6, 
                     seed=1000){
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
  
  # Make sure the collider is the cause of R-variable of its parent variable
  data_del <-gen_del(p, n, myDAG, mode, 
                     num_m=num_m, 
                     num_var=num_var, 
                     num_extra_e =num_extra_e )
  
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

gen_del <- function(p,n,myDAG,mode='mar',
                    num_m=6,
                    num_var=20, 
                    num_extra_e=3){
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
  prt_m = gen_m_prt(myDAG,mode,
                    num_m=num_m, 
                    num_var=num_var, 
                    num_extra_e =num_extra_e )
  
  
  # 3. Generate missing values according the the values of missingness indicators
  data_com_m_ref = generate_missing_values(p,n,myDAG, prt_m$ms, prt_m$prt_ms)
  data_m = data_com_m_ref$data_m
  data_ref = data_com_m_ref$data_ref
  data_complete = data_com_m_ref$data_complete
  
  return(list(data_complete=data_complete,
              data_m=data_m,
              data_ref=data_ref,
              parent_m_ind=prt_m$prt_ms,
              m_ind=prt_m$ms)) 
}


generate_missing_values <- function(p,n,myDAG,m_ind,parent_m_ind){
  # Give the parents of missingness indicators, and the missingness indcators
  # The missing values are generated whn the parent values are in the bottom XX percentage.
  data <- rmvnorm(n, mean=rep(0,p), sigma=trueCov(myDAG))
  data_m=data
  data_mcar=data
  for(i in c(1:length(m_ind))){
    # Choose lower "bottom_p" percentage of the values
    bottom_p <- runif(1, min = 0.1, max = 0.7)
    r <- data[,parent_m_ind[i]] < qnorm(bottom_p)
    data_m[r, m_ind[i]] = NA
    r_mcar<-sample(r)
    data_mcar[r_mcar, m_ind[i]] = NA
  }
  return(list(data_complete=data, 
              data_m=data_m,
              data_ref=data_mcar))
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
    cldr_prt[[count]]=a[m[,i]==1]
    count = count + 1 
  }
  cldr_prt
}

create_mnar_ind <- function(cldr,cldr_prt,num_var=20, num_extra_e=3, num_m = 6){
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
        if((pr !=cldr[i]) && (! (pr %in% ms))){
          ms[count] = pr
          prt_ms[count] = cldr[i]
          count = count + 1
        }
      }
    }
  }
  
  if(count > (num_extra_e+1)){
    ind_cld = sample(1:length(ms))
    ms = ms[ind_cld]
    prt_ms = prt_ms[ind_cld]
    count = (num_extra_e+1)
  }
  ind_rd = count
  left_ind = setdiff(1:num_var, ms)
  left_ind = sample(left_ind)
  for(i in 1:floor(length(left_ind)/2)){
    ms[count] = left_ind[i]
    prt_ms[count] = left_ind[floor(length(left_ind)/2)+i]
    count = count + 1
  }
  ms_ = ms[1:(ind_rd-1)]
  prt_ms_=prt_ms[1:(ind_rd-1)]
  ind = sample(ind_rd:length(ms))
  ms_[ind_rd:num_m] = ms[ind_rd:num_m]
  prt_ms_[ind_rd:num_m]=prt_ms[ind_rd:num_m]
  return(list(ms = ms_, prt_ms=prt_ms_))  
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
  
  if(count > (num_extra_e+1)){
    ind_cld = sample(1:length(ms))
    ms = ms[ind_cld]
    prt_ms = prt_ms[ind_cld]
    count = (num_extra_e+1)
  }
  ind_rd = count
  left_ind = setdiff(1:num_var, ms)
  left_ind = setdiff(left_ind, prt_ms)
  left_ind = sample(left_ind)
  for(i in 1:floor(length(left_ind)/2)){
    ms[count] = left_ind[i]
    prt_ms[count] = left_ind[floor(length(left_ind)/2)+i]
    count = count + 1
  }
  ms_ = ms[1:(ind_rd-1)]
  prt_ms_=prt_ms[1:(ind_rd-1)]
  ind = sample(ind_rd:length(ms))
  ms_[ind_rd:num_m] = ms[ind_rd:num_m]
  prt_ms_[ind_rd:num_m]=prt_ms[ind_rd:num_m]
  return(list(ms = ms_, prt_ms=prt_ms_))  
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

#****************** (Conditional) Independence Test ****************** 
binCItest.permc<- function(x, y, S, suffStat){  
  # print(c(x,y,S))
  if(!cond.PermC(x, y, S, suffStat)){return(binCItest_td(x,y,S,suffStat))}
  tryCatch(
    {
      ind_test <- c(x, y, S)
      ind_W <- get_prt_m_xys(c(x,y,S), suffStat)  # Get parents the {xyS} missingness indicators: prt_m
      ind_W <- setdiff(ind_W,ind_test)
      if(length(ind_W)==0){return(binCItest_td(x,y,S,suffStat))}
      ind_permc <- c(ind_test, ind_W)  
      
      ## Step 1: Get CPD given W = 0 and W = 1 
      data <- test_wise_deletion(ind_permc, suffStat$data)
      data <- data[,ind_permc] # Data are used for estimating Conditional Probability Distribution (CPD)
      num.td <- length(data[,1])
      num.W <- length(ind_W)
      comb.W <- bincombinations(length(ind_W))
      
      dm <- dim(comb.W)
      n.comb.W <- dm[1]
      
      p.joint.w <- list()
      for(i in 1:n.comb.W){
        ind.mask <- comb.W[i, ] == data[,(length(ind_test)+1):length(ind_permc)]
        ind.mask <- matrix(data = ind.mask, nrow = length(data[,1]))
        ind_wi <- apply(ind.mask, 1, function(x) sum(x)==num.W)  # every element is the same with W pattern
        dat.wi <- data[ind_wi, 1:length(ind_test)]
        p.xy.wi <- apply(dat.wi,2, function(x) sum(x)/length(x))
        cor.xy.wi <- cor(dat.wi)
        p.joint.w[[i]] <- ObtainMultBinaryDist(corr = cor.xy.wi, marg.probs = p.xy.wi )
      }
      
      ## Step 2: Shuffle W 
      data.W <- test_wise_deletion(ind_W, suffStat$data)
      nrow.td.W <- length(data.W[,1])
      data.W <- data.W[,ind_W]
      if(length(ind_W)==1){
        data.W <- matrix(data = data.W, nrow = nrow.td.W,byrow = TRUE)  
      }
      
      ind <- 1:nrow.td.W
      ind.perm.w <- sample(ind)[1:num.td]
      data.W.perm <- data.W[ind.perm.w,]
      if(length(ind_W)==1){data.W.perm <- matrix(data = data.W.perm, nrow = num.td,byrow = TRUE)}
      
      ## Step 3: Generate virtual data of x,y, and S
      for(i in 1:n.comb.W){
        ind.mask <- comb.W[i, ] == data.W.perm
        if(length(ind_W)==1){ind.mask <- matrix(data = ind.mask, nrow = length(data.W.perm[,1]),byrow = TRUE)}
        
        ind_wi <- apply(ind.mask, 1, function(x) sum(x)==num.W)  # every element is the same with W pattern
        
        if(i==1){
          data.vir <- RMultBinary(n = sum(ind_wi), mult.bin.dist = p.joint.w[[i]])$binary.sequences  
        }else{
          data.vir <- rbind(data.vir, RMultBinary(n = sum(ind_wi), mult.bin.dist = p.joint.w[[i]])$binary.sequences  )
        }
      }
      
      ## Step4: Test with the virtual data set
      if(length(S) > 0){binCItest(1,2,c(3:length(ind_test)), list(dm = data.vir, adaptDF = TRUE))}
      else{binCItest(1,2,c(), list(dm = data.vir , adaptDF = TRUE))}
    },
    error=function(cond) {
      message(cond)
      # Choose a return value in case of error
      return(binCItest_td(x,y,S,suffStat))
    })
    
}  

binCItest.drw <- function(x, y, S, suffStat){
  weights_ = compute_weights(x, y, S, suffStat)
  test_ind = c(x,y,S)
  del_res = test_wise_deletion_w(test_ind,suffStat$data, weights_)
  
  data = del_res$data
  weights = del_res$weights
  if(length(S)>0){
    binCItest_w(1, 2, 3:length(test_ind), weights,list(dm = data[,test_ind], adaptDF = FALSE))
  }
  else{
    binCItest_w(1, 2, c(),weights, list(dm = data[,test_ind], adaptDF = FALSE))
  }
  
}

binCItest_td_ref<- function(x, y, S, suffStat){
  test_ind = c(x,y,S)
  data = test_wise_deletion(test_ind,suffStat$data)
  sample_size <<- c(sample_size,length(data[,1]))
  if(length(S)>0){
    binCItest(1, 2, 3:length(test_ind), list(dm = data[,test_ind], adaptDF = FALSE))  
  }
  else{
    binCItest(1, 2, c(), list(dm = data[,test_ind], adaptDF = FALSE))
  }
}

binCItest_td <- function(x, y, S, suffStat){
  test_ind = c(x,y,S)
  data = test_wise_deletion(test_ind,suffStat$data)
  # sample_size <<- c(sample_size,length(data[,1]))
  if(length(S)>0){
    binCItest(1, 2, 3:length(test_ind), list(dm = data[,test_ind], adaptDF = FALSE))  
  }
  else{
    binCItest(1, 2, c(), list(dm = data[,test_ind], adaptDF = FALSE))
  }
  
}

f_R<-function(r,suffStat){
  n.sample = dim(suffStat$data)[1]
  if(length(r) > 1){
    R = suffStat$data[,r]
  }
  else{R = matrix(data = suffStat$data[,r], nrow = n.sample) # BUG: here , only get one column
  }
  R = is.na(R) * 1
  return(sum(rowSums(R)==0)/n.sample)
}

f_weights<-function(ri, pa, pval, suffStat){
  data = suffStat$data[,pa]
  ri = is.na(suffStat$data[,ri]) * 1
  dat_r = cbind(ri,data)
  n.var = dim(dat_r)[2]
  dat_r = test_wise_deletion(1:n.var, dat_r)
  return(sum(dat_r[,1] == 0 & dat_r[,2:n.var]==pval)/sum( dat_r[,2:n.var]==pval ))
}

get_rw_pair<-function(x,y,S,ind_W,suffStat){
  rw <- list(r=c(),w=list())
  
  count = 1
  for(ri in c(x,y,S,ind_W)){
    wi <- get_prt_m_xys(ri, suffStat)
    if(length(wi)!=0){
      rw$r[count] = ri
      rw$w[[count]] = wi
      count = count + 1
    }
  }
  rw
}

weights_check_table<-function(rw,suffStat,indW, comb.W){
  weights_tab = c()
  pR <- f_R(rw$r,suffStat)
  n.comb.W <- dim(comb.W)[1]
  for(i in 1:n.comb.W){
    weighti = pR
    for(count in 1:length(rw$r)){
      ri = rw$r[count]
      wi =rw$w[[count]]
      pos_wi = (wi == indW) %*% 1:length(indW)
      val_wi = comb.W[i, pos_wi]
      weighti = weighti / f_weights(ri, wi, val_wi, suffStat)
    }
    weights_tab = c(weights_tab, weighti)
  }
  weights_tab
}

compute_weights<- function(x, y, S, suffStat){
  data = suffStat$data
  n.sample = dim(data)[1]
  weights = rep(1, n.sample)
  
  # Detection of parents of missingness indicators
  ind_test <- c(x, y, S)
  ind_W <- unique(get_prt_m_xys(c(x,y,S), suffStat))  # Get parents the {xyS} missingness indicators
  if(length(ind_W)==0){return(weights)}
  
  checkW <- ind_W
  while(length(pa_W <- get_prt_m_xys(checkW, suffStat)) > 0){
    ind_W <- c(ind_W, pa_W) # Get parents the W missingness indicators
    checkW = pa_W
  } 
  ind_W <- unique(ind_W)
  # Get ri and corresponding wi
  rw = get_rw_pair(x,y,S,ind_W,suffStat)
  # Get the weights check table weights <==> W
  num.W <- length(ind_W)
  comb.W <- bincombinations(length(ind_W))
  weights_tab = weights_check_table(rw,suffStat,ind_W,comb.W)
  # apply the weights check table to assigning weights for each data point
  if(length(ind_W)>1 ){
    W = data[,ind_W]
  }
  else{
    W = matrix(data = data[,ind_W], nrow = n.sample)  
  }
  
  for(i in 1: dim(comb.W)[1]){
    ind = rep(TRUE, n.sample)
    for(j in 1:dim(comb.W)[2]){
      ind = (W[,j] == comb.W[i,j]) & ind
    }
    weights[ind] =  weights_tab[i] 
  }
  return(weights)
}

binCItest_w <- function (x, y, S, weights, suffStat)
{
  if (is.data.frame(dm <- suffStat$dm))
    dm <- data.matrix(dm)
  else stopifnot(is.matrix(dm))
  adaptDF <- suffStat$adaptDF
  gSquareBin_w(x = x, y = y, S = S, dm = dm, weights=weights, adaptDF = adaptDF,
               verbose = FALSE)
}

gSquareBin_w <- function (x, y, S, dm, weights, adaptDF = FALSE, n.min = 10 * df, verbose = FALSE)
{
  stopifnot((n <- nrow(dm)) >= 1)
  if (!all(as.logical(dm) == dm))
    stop("'dm' must be binary, i.e. with values in {0,1}")
  if (verbose)
    cat("Edge ", x, "--", y, " with subset S =", S, "\n")
  lenS <- length(S)
  df <- 2^lenS
  if (n < n.min) {
    warning(gettextf("n=%d is too small (n < n.min = %d ) for G^2 test (=> treated as independence)",
                     n, n.min), domain = NA)
    return(1)
  }
  d.x1 <- dm[, x] + 1L
  d.y1 <- dm[, y] + 1L
  if (lenS <= 5) {
    n12 <- 1:2
    switch(lenS + 1L, {
      nijk <- array(0L, c(2, 2))
      for (i in n12) {
        d.x.i <- d.x1 == i
        for (j in n12) nijk[i, j] <- sum(weights[d.x.i & d.y1 ==j])
      }
      t.log <- n * (nijk/tcrossprod(rowSums(nijk), colSums(nijk)))
    }, {
      dmS.1 <- dm[, S] + 1L
      nijk <- array(0L, c(2, 2, 2))
      for (i in n12) {
        d.x.i <- d.x1 == i
        for (j in n12) {
          d.x.i.y.j <- d.x.i & d.y1 == j
          for (k in n12) nijk[i, j, k] <- sum(weights[d.x.i.y.j &
                                                        dmS.1 == k])
        }
      }
      alt <- c(x, y, S)
      c <- which(alt == S)
      nik <- apply(nijk, c, rowSums)
      njk <- apply(nijk, c, colSums)
      nk <- colSums(njk)
      t.log <- array(0, c(2, 2, 2))
      if (c == 3) {
        for (k in n12) t.log[, , k] <- nijk[, , k] *
            (nk[k]/tcrossprod(nik[, k], njk[, k]))
      } else if (c == 1) {
        for (k in n12) t.log[k, , ] <- nijk[k, , ] *
            (nk[k]/tcrossprod(nik[, k], njk[, k]))
      } else {
        for (k in n12) t.log[, k, ] <- nijk[, k, ] *
            (nk[k]/tcrossprod(nik[, k], njk[, k]))
      }
    }, {
      dmS1.1 <- dm[, S[1]] + 1L
      dmS2.1 <- dm[, S[2]] + 1L
      nijk2 <- array(0L, c(2, 2, 2, 2))
      for (i in n12) {
        d.x.i <- d.x1 == i
        for (j in n12) {
          d.x.i.y.j <- d.x.i & d.y1 == j
          for (k in n12) {
            d.x.y.S1 <- d.x.i.y.j & dmS1.1 == k
            for (l in n12) nijk2[i, j, k, l] <- sum(weights[d.x.y.S1 &
                                                              dmS2.1 == l])
          }
        }
      }
      nijk <- array(0L, c(2, 2, 4))
      for (i in n12) {
        for (j in n12) nijk[, , 2 * (i - 1) + j] <- nijk2[,
                                                          , i, j]
      }
      nik <- apply(nijk, 3, rowSums)
      njk <- apply(nijk, 3, colSums)
      nk <- colSums(njk)
      t.log <- array(0, c(2, 2, 4))
      for (k in 1:4) t.log[, , k] <- nijk[, , k] * (nk[k]/tcrossprod(nik[,
                                                                         k], njk[, k]))
    }, {
      dmS1.1 <- dm[, S[1]] + 1L
      dmS2.1 <- dm[, S[2]] + 1L
      dmS3.1 <- dm[, S[3]] + 1L
      nijk <- array(0L, c(2, 2, 8))
      for (i1 in n12) {
        d.x.i <- d.x1 == i1
        for (i2 in n12) {
          d.x.y.i12 <- d.x.i & d.y1 == i2
          for (i3 in n12) {
            d.x.y.S1 <- d.x.y.i12 & dmS1.1 == i3
            for (i4 in n12) {
              d.xy.S1S2 <- d.x.y.S1 & dmS2.1 == i4
              for (i5 in n12) nijk[i1, i2, 4 * (i3 -
                                                  1) + 2 * (i4 - 1) + i5] <- sum(weights[d.xy.S1S2 &
                                                                                           dmS3.1 == i5])
            }
          }
        }
      }
      nik <- apply(nijk, 3, rowSums)
      njk <- apply(nijk, 3, colSums)
      nk <- colSums(njk)
      t.log <- array(0, c(2, 2, 8))
      for (k in 1:8) t.log[, , k] <- nijk[, , k] * (nk[k]/tcrossprod(nik[,
                                                                         k], njk[, k]))
    }, {
      dmS1.1 <- dm[, S[1]] + 1L
      dmS2.1 <- dm[, S[2]] + 1L
      dmS3.1 <- dm[, S[3]] + 1L
      dmS4.1 <- dm[, S[4]] + 1L
      nijk <- array(0L, c(2, 2, 16))
      for (i1 in n12) {
        d.x.i <- d.x1 == i1
        for (i2 in n12) {
          d.x.y.i12 <- d.x.i & d.y1 == i2
          for (i3 in n12) {
            d.x.y.S1 <- d.x.y.i12 & dmS1.1 == i3
            for (i4 in n12) {
              d.xy.S1S2 <- d.x.y.S1 & dmS2.1 == i4
              for (i5 in n12) {
                d.xy.S1S2S3 <- d.xy.S1S2 & dmS3.1 ==
                  i5
                for (i6 in n12) nijk[i1, i2, 8 * (i3 -
                                                    1) + 4 * (i4 - 1) + 2 * (i5 - 1) +
                                       i6] <- sum(weights[d.xy.S1S2S3 & dmS4.1 ==
                                                            i6])
              }
            }
          }
        }
      }
      nik <- apply(nijk, 3, rowSums)
      njk <- apply(nijk, 3, colSums)
      nk <- colSums(njk)
      t.log <- array(0, c(2, 2, 16))
      for (k in 1:16) t.log[, , k] <- nijk[, , k] * (nk[k]/tcrossprod(nik[,
                                                                          k], njk[, k]))
    }, {
      dmS1.1 <- dm[, S[1]] + 1L
      dmS2.1 <- dm[, S[2]] + 1L
      dmS3.1 <- dm[, S[3]] + 1L
      dmS4.1 <- dm[, S[4]] + 1L
      dmS5.1 <- dm[, S[5]] + 1L
      nijk <- array(0L, c(2, 2, 32))
      for (i1 in n12) {
        d.x.i <- d.x1 == i1
        for (i2 in n12) {
          d.x.y.i12 <- d.x.i & d.y1 == i2
          for (i3 in n12) {
            d.x.y.S1 <- d.x.y.i12 & dmS1.1 == i3
            for (i4 in n12) {
              d.xy.S1S2 <- d.x.y.S1 & dmS2.1 == i4
              for (i5 in n12) {
                d.xy.S1S2S3 <- d.xy.S1S2 & dmS3.1 ==
                  i5
                for (i6 in n12) {
                  d.xy.S1S2S3S4 <- d.xy.S1S2S3 & dmS4.1 ==
                    i6
                  for (i7 in n12) nijk[i1, i2, 16 *
                                         (i3 - 1) + 8 * (i4 - 1) + 4 * (i5 -
                                                                          1) + 2 * (i6 - 1) + i7] <- sum(weights[d.xy.S1S2S3S4 &
                                                                                                                   dmS5.1 == i7])
                }
              }
            }
          }
        }
      }
      nik <- apply(nijk, 3, rowSums)
      njk <- apply(nijk, 3, colSums)
      nk <- colSums(njk)
      t.log <- array(0, c(2, 2, 32))
      for (k in 1:32) t.log[, , k] <- nijk[, , k] * (nk[k]/tcrossprod(nik[,
                                                                          k], njk[, k]))
    })
  }
  else {
    nijk <- array(0L, c(2, 2, 1))
    i <- d.x1[1]
    j <- d.y1[1]
    k <- NULL
    lapply(as.list(S), function(x) {
      k <<- cbind(k, dm[, x] + 1)
      NULL
    })
    parents.count <- 1L
    parents.val <- t(k[1, ])
    nijk[i, j, parents.count] <- 1L
    for (it.sample in 2:n) {
      new.p <- TRUE
      i <- d.x1[it.sample]
      j <- d.y1[it.sample]
      t.comp <- t(parents.val[1:parents.count, ]) == k[it.sample,
                                                       ]
      dim(t.comp) <- c(lenS, parents.count)
      for (it.parents in 1:parents.count) {
        if (all(t.comp[, it.parents])) {
          nijk[i, j, it.parents] <- nijk[i, j, it.parents] +
            1L
          new.p <- FALSE
          break
        }
      }
      if (new.p) {
        parents.count <- parents.count + 1L
        if (verbose >= 2)
          cat(sprintf(" adding new parents (count = %d) at sample %d\n",
                      parents.count, it.sample))
        parents.val <- rbind(parents.val, k[it.sample,
                                            ])
        nijk <- abind(nijk, array(0, c(2, 2, 1)))
        nijk[i, j, parents.count] <- 1L
      }
    }
    nik <- apply(nijk, 3, rowSums)
    njk <- apply(nijk, 3, colSums)
    nk <- colSums(njk)
    t.log <- array(0, c(2, 2, parents.count))
    for (k in 1:parents.count) t.log[, , k] <- nijk[, ,
                                                    k] * (nk[k]/tcrossprod(nik[, k], njk[, k]))
  }
  G2 <- sum(2 * nijk * log(t.log), na.rm = TRUE)
  if (adaptDF && lenS > 0) {
    zero.counts <- sum(nijk == 0L) + 4 * (2^lenS - dim(nijk)[3])
    ndf <- max(1, df - zero.counts)
    if (verbose)
      cat("adaptDF: (df=", df, ", zero.counts=", zero.counts,
          ") ==> new df = ", ndf, "\n", sep = "")
    df <- ndf
  }
  pchisq(G2, df, lower.tail = FALSE)
}

gSquareBin.weighted <- function (x, y, S, W, dm, adaptDF = FALSE, n.min = 10 * df, verbose = FALSE) 
{
  stopifnot((n <- nrow(dm)) >= 1)
  if (!all(as.logical(dm) == dm)) 
    stop("'dm' must be binary, i.e. with values in {0,1}")
  if (verbose) 
    cat("Edge ", x, "--", y, " with subset S =", S, "\n")
  lenS <- length(S)
  df <- 2^lenS
  if (n < n.min) {
    warning(gettextf("n=%d is too small (n < n.min = %d ) for G^2 test (=> treated as independence)", 
                     n, n.min), domain = NA)
    return(1)
  }
  d.x1 <- dm[, x] + 1L
  d.y1 <- dm[, y] + 1L
  if (lenS <= 5) {
    n12 <- 1:2
    switch(lenS + 1L, {
      nijk <- array(0L, c(2, 2))
      for (i in n12) {
        d.x.i <- d.x1 == i
        for (j in n12) nijk[i, j] <- sum(d.x.i & d.y1 == 
                                           j)
      }
      t.log <- n * (nijk/tcrossprod(rowSums(nijk), colSums(nijk)))
    }, {
      dmS.1 <- dm[, S] + 1L
      nijk <- array(0L, c(2, 2, 2))
      for (i in n12) {
        d.x.i <- d.x1 == i
        for (j in n12) {
          d.x.i.y.j <- d.x.i & d.y1 == j
          for (k in n12) nijk[i, j, k] <- sum(d.x.i.y.j & 
                                                dmS.1 == k)
        }
      }
      alt <- c(x, y, S)
      c <- which(alt == S)
      nik <- apply(nijk, c, rowSums)
      njk <- apply(nijk, c, colSums)
      nk <- colSums(njk)
      t.log <- array(0, c(2, 2, 2))
      if (c == 3) {
        for (k in n12) t.log[, , k] <- nijk[, , k] * 
            (nk[k]/tcrossprod(nik[, k], njk[, k]))
      } else if (c == 1) {
        for (k in n12) t.log[k, , ] <- nijk[k, , ] * 
            (nk[k]/tcrossprod(nik[, k], njk[, k]))
      } else {
        for (k in n12) t.log[, k, ] <- nijk[, k, ] * 
            (nk[k]/tcrossprod(nik[, k], njk[, k]))
      }
    }, {
      dmS1.1 <- dm[, S[1]] + 1L
      dmS2.1 <- dm[, S[2]] + 1L
      nijk2 <- array(0L, c(2, 2, 2, 2))
      for (i in n12) {
        d.x.i <- d.x1 == i
        for (j in n12) {
          d.x.i.y.j <- d.x.i & d.y1 == j
          for (k in n12) {
            d.x.y.S1 <- d.x.i.y.j & dmS1.1 == k
            for (l in n12) nijk2[i, j, k, l] <- sum(d.x.y.S1 & 
                                                      dmS2.1 == l)
          }
        }
      }
      nijk <- array(0L, c(2, 2, 4))
      for (i in n12) {
        for (j in n12) nijk[, , 2 * (i - 1) + j] <- nijk2[, 
                                                          , i, j]
      }
      nik <- apply(nijk, 3, rowSums)
      njk <- apply(nijk, 3, colSums)
      nk <- colSums(njk)
      t.log <- array(0, c(2, 2, 4))
      for (k in 1:4) t.log[, , k] <- nijk[, , k] * (nk[k]/tcrossprod(nik[, 
                                                                         k], njk[, k]))
    }, {
      dmS1.1 <- dm[, S[1]] + 1L
      dmS2.1 <- dm[, S[2]] + 1L
      dmS3.1 <- dm[, S[3]] + 1L
      nijk <- array(0L, c(2, 2, 8))
      for (i1 in n12) {
        d.x.i <- d.x1 == i1
        for (i2 in n12) {
          d.x.y.i12 <- d.x.i & d.y1 == i2
          for (i3 in n12) {
            d.x.y.S1 <- d.x.y.i12 & dmS1.1 == i3
            for (i4 in n12) {
              d.xy.S1S2 <- d.x.y.S1 & dmS2.1 == i4
              for (i5 in n12) nijk[i1, i2, 4 * (i3 - 
                                                  1) + 2 * (i4 - 1) + i5] <- sum(d.xy.S1S2 & 
                                                                                   dmS3.1 == i5)
            }
          }
        }
      }
      nik <- apply(nijk, 3, rowSums)
      njk <- apply(nijk, 3, colSums)
      nk <- colSums(njk)
      t.log <- array(0, c(2, 2, 8))
      for (k in 1:8) t.log[, , k] <- nijk[, , k] * (nk[k]/tcrossprod(nik[, 
                                                                         k], njk[, k]))
    }, {
      dmS1.1 <- dm[, S[1]] + 1L
      dmS2.1 <- dm[, S[2]] + 1L
      dmS3.1 <- dm[, S[3]] + 1L
      dmS4.1 <- dm[, S[4]] + 1L
      nijk <- array(0L, c(2, 2, 16))
      for (i1 in n12) {
        d.x.i <- d.x1 == i1
        for (i2 in n12) {
          d.x.y.i12 <- d.x.i & d.y1 == i2
          for (i3 in n12) {
            d.x.y.S1 <- d.x.y.i12 & dmS1.1 == i3
            for (i4 in n12) {
              d.xy.S1S2 <- d.x.y.S1 & dmS2.1 == i4
              for (i5 in n12) {
                d.xy.S1S2S3 <- d.xy.S1S2 & dmS3.1 == 
                  i5
                for (i6 in n12) nijk[i1, i2, 8 * (i3 - 
                                                    1) + 4 * (i4 - 1) + 2 * (i5 - 1) + 
                                       i6] <- sum(d.xy.S1S2S3 & dmS4.1 == 
                                                    i6)
              }
            }
          }
        }
      }
      nik <- apply(nijk, 3, rowSums)
      njk <- apply(nijk, 3, colSums)
      nk <- colSums(njk)
      t.log <- array(0, c(2, 2, 16))
      for (k in 1:16) t.log[, , k] <- nijk[, , k] * (nk[k]/tcrossprod(nik[, 
                                                                          k], njk[, k]))
    }, {
      dmS1.1 <- dm[, S[1]] + 1L
      dmS2.1 <- dm[, S[2]] + 1L
      dmS3.1 <- dm[, S[3]] + 1L
      dmS4.1 <- dm[, S[4]] + 1L
      dmS5.1 <- dm[, S[5]] + 1L
      nijk <- array(0L, c(2, 2, 32))
      for (i1 in n12) {
        d.x.i <- d.x1 == i1
        for (i2 in n12) {
          d.x.y.i12 <- d.x.i & d.y1 == i2
          for (i3 in n12) {
            d.x.y.S1 <- d.x.y.i12 & dmS1.1 == i3
            for (i4 in n12) {
              d.xy.S1S2 <- d.x.y.S1 & dmS2.1 == i4
              for (i5 in n12) {
                d.xy.S1S2S3 <- d.xy.S1S2 & dmS3.1 == 
                  i5
                for (i6 in n12) {
                  d.xy.S1S2S3S4 <- d.xy.S1S2S3 & dmS4.1 == 
                    i6
                  for (i7 in n12) nijk[i1, i2, 16 * 
                                         (i3 - 1) + 8 * (i4 - 1) + 4 * (i5 - 
                                                                          1) + 2 * (i6 - 1) + i7] <- sum(d.xy.S1S2S3S4 & 
                                                                                                           dmS5.1 == i7)
                }
              }
            }
          }
        }
      }
      nik <- apply(nijk, 3, rowSums)
      njk <- apply(nijk, 3, colSums)
      nk <- colSums(njk)
      t.log <- array(0, c(2, 2, 32))
      for (k in 1:32) t.log[, , k] <- nijk[, , k] * (nk[k]/tcrossprod(nik[, 
                                                                          k], njk[, k]))
    })
  }
  else {
    nijk <- array(0L, c(2, 2, 1))
    i <- d.x1[1]
    j <- d.y1[1]
    k <- NULL
    lapply(as.list(S), function(x) {
      k <<- cbind(k, dm[, x] + 1)
      NULL
    })
    parents.count <- 1L
    parents.val <- t(k[1, ])
    nijk[i, j, parents.count] <- 1L
    for (it.sample in 2:n) {
      new.p <- TRUE
      i <- d.x1[it.sample]
      j <- d.y1[it.sample]
      t.comp <- t(parents.val[1:parents.count, ]) == k[it.sample, 
                                                       ]
      dim(t.comp) <- c(lenS, parents.count)
      for (it.parents in 1:parents.count) {
        if (all(t.comp[, it.parents])) {
          nijk[i, j, it.parents] <- nijk[i, j, it.parents] + 
            1L
          new.p <- FALSE
          break
        }
      }
      if (new.p) {
        parents.count <- parents.count + 1L
        if (verbose >= 2) 
          cat(sprintf(" adding new parents (count = %d) at sample %d\n", 
                      parents.count, it.sample))
        parents.val <- rbind(parents.val, k[it.sample, 
                                            ])
        nijk <- abind(nijk, array(0, c(2, 2, 1)))
        nijk[i, j, parents.count] <- 1L
      }
    }
    nik <- apply(nijk, 3, rowSums)
    njk <- apply(nijk, 3, colSums)
    nk <- colSums(njk)
    t.log <- array(0, c(2, 2, parents.count))
    for (k in 1:parents.count) t.log[, , k] <- nijk[, , 
                                                    k] * (nk[k]/tcrossprod(nik[, k], njk[, k]))
  }
  G2 <- sum(2 * nijk * log(t.log), na.rm = TRUE)
  if (adaptDF && lenS > 0) {
    zero.counts <- sum(nijk == 0L) + 4 * (2^lenS - dim(nijk)[3])
    ndf <- max(1, df - zero.counts)
    if (verbose) 
      cat("adaptDF: (df=", df, ", zero.counts=", zero.counts, 
          ") ==> new df = ", ndf, "\n", sep = "")
    df <- ndf
  }
  pchisq(G2, df, lower.tail = FALSE)
}

gaussCItest_td <- function(x, y, S, suffStat) {
  ## Conditional independence test between continuous variables with deletion methods
  ## test P(x,y|S)
  ## suffStat: the class contains the dataset and other necessary variables
  ##    suffStat$data: the dataset
  ##--------------
  ## Return: the p-value of the test 
  
  data = test_wise_deletion(c(x,y,S),suffStat$data)
  suffStat$data = data
  suffStat$C = cor(data)
  suffStat$n = length(data[,1])
  gaussCItest(x, y, S, suffStat)
}

gaussCItest.drw <- function(x, y, S, suffStat) {
  ## Conditional independence test between continuous variables with deletion methods
  ## test P(x,y|S)
  ## suffStat: the class contains the dataset and other necessary variables
  ##    suffStat$data: the dataset
  ##--------------
  ## Return: the p-value of the test 
  # Note that "tw_data" and "weights" should have the same order/index
  
  ind_W <- unique(get_prt_m_xys(c(x,y,S), suffStat))  # Get parents the {xyS} missingness indicators
  if(length(ind_W)==0){return(gaussCItest_td(x,y,S,suffStat))}
  
  checkW <- ind_W
  while(length(pa_W <- get_prt_m_xys(checkW, suffStat)) > 0){
    ind_W <- c(ind_W, pa_W) # Get parents the W missingness indicators
    checkW = pa_W
  } 
  ind_W <- unique(ind_W)
  
  corr_ind = c(x,y,S,ind_W)
  tw_data = test_wise_deletion(corr_ind, suffStat$data)
  weights = compute.weights.continuous(corr_ind, suffStat)
 
  suffStat$C = wtd.cors(tw_data, tw_data, weights)
  suffStat$n = length(weights)
  
  if(length(S)>0){
    gaussCItest(1, 2, 3:length(c(x,y,S)),suffStat)
  }
  else{
    gaussCItest(1, 2, c(),suffStat)
  }
}

compute.weights.continuous<-function(corr_ind, suffStat){
  ind.twdel = indx_test_wise_deletion(corr_ind,suffStat$data)
  weights = rep(1,length(ind.twdel))  # length of test-wise deleted data
  ind_r_xys = get_ind_r_xys(corr_ind, suffStat)
  
  for(ind_ri in ind_r_xys){
    ind_pa =get_prt_m_xys(c(ind_ri), suffStat) 
    # Return the single parent of a missingness indicator (an element, i.e., list[[i]], not a list, i.e., list[i])
    pa = suffStat$data[,ind_pa]
    beta = k.weights(pa[ind.twdel], pa[!is.na(pa)])
    weights = weights * beta
  }
  weights
}

k.weights <- function (x, x_target){
  # kernel-based density ratio estimate
  setVariable(matlab, X = x)
  setVariable(matlab, X_target=x_target)
  if(length(x)<800){
    coeff = 0.5 
  }else if(length(x) > 1500){
    coeff = 0.25
  }else{
    coeff = 0.4
  }
  set.width = paste0("width =",coeff,"* std(X);", sep = "")
  
  evaluate(matlab, set.width)
  evaluate(matlab, "[beta_cs EXITFLAG_cs] = betaKMM_improved(X, X_target, width, 0, 0);")
  beta_cs = getVariable(matlab, "beta_cs")
  beta = beta_cs$beta.cs
}

PermCCItest <- function(x, y, S, suffStat){
  ## The Z <- XY; Rz <- XY is not included in the test 
  ## Step 1: Learning generaive model for {X, Y, S} to impute X, Y, and S
  if(!cond.PermC(x, y, S, suffStat)){return(gaussCItest_td(x,y,S,suffStat))}
  ind_W = get_prt_m_xys(c(x,y,S), suffStat)  # Get parents the {xyS} missingness indicators: prt_m
  ind_permc <- c(x, y, S, ind_W)
  ind_test <- c(x, y, S)
  data <- test_wise_deletion(ind_permc, suffStat$data)
  data <- data[,ind_permc]
  
  xnam <- paste0("data[,", (length(ind_test)+1):length(ind_permc),"]")
  lr <- list()
  res <- list()
  for( i in 1:length(ind_test)){  # NOTE: the order always the "1(x) CI 2(y) given all the others (S)"  The order of xyS and the order of PermC test 
    fmla <- as.formula(paste("data[,", i,"] ~ 0 + ", paste(xnam, collapse= "+"),""))
    lr[[i]] <- lm(formula = fmla, data = data.frame(data))  # the ith block of list is ...
    res[[i]] <- residuals(lr[[i]])  # residuals
  }
  
  ## Step 2: Shuffle the source "W" -- parents of the missingness indicators
  #         a) Remove the missing entries of dat[, W]
  #         b) Row-based/ sample based Shuffle the index of W data points
  #         c) the same number of test-wise deletion CI Test 
  data_W_p = perm(ind_W, suffStat$data)
  for(i in (length(ind_test)+1):length(ind_permc)){
    data[, i] = data_W_p[1:length(data[,1]),i-length(ind_test)]  
  }
  
  
  ## Step 3: Generate the virtual data follows the full data distribution P(X, Y, S)
  vir <- list()
  for(i in 1:length(ind_test)){  # NOTE: the order always the "1(x) CI 2(y) given all the others (S)"  The order of xyS and the order of PermC test 
    vir[[i]] = predict(lr[[i]],list(data)) + res[[i]] 
  }
  
  data_perm = data.frame(vir[[1]])
  for(i in 1:length(ind_test)){ # NOTE: the order always the "1(x) CI 2(y) given all the others (S)"  The order of xyS and the order of PermC test 
    data_perm[,i]=vir[[i]]
  }
  suffStat_perm = list(C=cor(data_perm), n=length(data_perm[,1]))
  if(length(ind_test)>2){
    gaussCItest(1, 2, 3:length(ind_test), suffStat_perm)  
  }
  else{
    gaussCItest(1, 2,c(), suffStat_perm)
  }
}

DRWCItest <- function(x, y, S, suffStat){
  
  # Determine whether the current case satisfies the condition of correction
  if(!cond.PermC(x, y, S, suffStat)){return(gaussCItest_td(x,y,S,suffStat))}
  
  # Logistic regression for each missingness indicator (rx, ry, rsi): ri ~ Pa(ri).
  # Note that data for training the logistic regression are different with the data used for computing weights.
  ind_r_xys <- get_ind_r_xys(c(x, y, S), suffStat)
  ind <- get_ind_weights(ind_r_xys,  suffStat)  # ind is a list of logical variable
  weights = sum(ind)/length(ind)  # "ri = 0" means missing value
  
  for(ind_ri in ind_r_xys){
    logidata <- get_logidata(ind_ri,  suffStat)
    fmla <- get_logi_formula(length(logidata))
    logi.fit <- glm(formula = as.formula(fmla), data = logidata, family = binomial)
    
    prt_i <- get_prt_i(ind_ri, suffStat)
    
    logidata = data.frame(r = 1:length(sum(ind)), prt = suffStat$data[ind,prt_i])
    logi.prob = predict(logi.fit, logidata, type = "response")  
    # TODO: check what if the name of the column is different with the formula
    # it should have the same name with the formula
    # it follow the formula slicing way when input the same name with formula
    weights = weights * (1/logi.prob)
  }
  
  # Weighted the correlation
  data_td <- suffStat$data[ind,]
  suffStat_drw  = list(C=wtd.cors(data_td,data_td,weights), n=length(weights))
  gaussCItest(x, y, S, suffStat_drw)
  
}  
# Logistic regression is not general enough for estimate the ratio, thus waiting for density ratio estimate

iscorr<- function(x, y, S, suffStat){
  return(TRUE)
}

get_logi_formula <- function(len){
  xnam <- paste0("logidata[,", 2:len,"]","")
  # as.formula(paste("logidata[,1] ~ ", paste(xnam, collapse= "+")))
  paste("logidata[,1] ~ ", paste(xnam, collapse= "+"))
}

get_ind_weights <- function(ind_r_xys,  suffStat){
  var_ind <- c(ind_r_xys)
  for(ind_ri in ind_r_xys){
    var_ind <- union(var_ind,get_prt_i(ind_ri, suffStat))
  }
  data = suffStat$data
  
  not_del_ind = c(rep(TRUE,length(data[,1])))
  for (var in var_ind){
    if(anyNA(data[,var])){
      not_del_ind = not_del_ind & !is.na(data[,var])
    }
  } 
  return(not_del_ind)
}

get_prt_i<- function(ind_ri, suffStat){
  prt_cell <- suffStat$prt_m['prt'][suffStat$prt_m['m'] == ind_ri]
  prt_cell[[1]]
}
  
get_logidata <- function(ind_ri,  suffStat){
  prt_i <- get_prt_i(ind_ri, suffStat)
  ri <- as.integer(!is.na(suffStat$data[,ind_ri]))
  logidata <- data.frame(ri, suffStat$data[,prt_i])
  test_wise_deletion(1:length(logidata[1,]), logidata)  
}

get_ind_r_xys <- function(ind, suffStat){
  # return index of missingness indicators of x, y, S
  ind_r_xys <- c()
  for(i in ind){
    if(sum(is.na(suffStat$data[,i]))>0){
      ind_r_xys = c(ind_r_xys,i) 
      }
  }
  return(ind_r_xys)
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

indx_test_wise_deletion <-function(var_ind, data){
  ## Delete the rows of given variables (var_ind) if there is a missing value in a row
  ## var_ind: variables in the current conditonal independence test
  ## data: the whole data set 
  ##--------------------------------
  ## Return: the index of the deleted dataset with the variables in the CI test 
  
  not_del_ind = c(rep(TRUE,length(data[,1])))
  for (var in var_ind){
    if(anyNA(data[,var])){
      not_del_ind = not_del_ind & !is.na(data[,var])
    }
  }  
  ind = 1:length(not_del_ind)
  return(ind[not_del_ind])
}

test_wise_deletion_w <-function(var_ind, data,weights){
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
  not_del_ind = not_del_ind & !is.na(weights)
  return(list(data= data[not_del_ind,],weights = weights[not_del_ind]))
}

get_prt_m_xys<-function(ind, suffStat){
  w = c()
  for(i in ind){
    if(is.in_prt_m(i, suffStat$prt_m)){
      prt_i <- get_prt_i(i, suffStat)
      w = c(w, prt_i)  
    }
  }
  w
}

perm <- function(W, data){
  data <- test_wise_deletion(W,data)
  len = length(data[,1])
  ind_p <- sample(1:len)
  data = data[ind_p,]
  data_permw <- data[,W]
  data.frame(data_permw)
}

is.in_prt_m<-function(i, prt_m){
  sum(prt_m['m'] == i) > 0
}

cond.PermC<-function(x, y, S, suffStat){
  ind <- c(x,y,S)
  cond <- FALSE
  if("skel" %in% names(suffStat)){
    if(length(intersect(ind, suffStat$prt_m$m)) > 0){  # 1) xyS have missingness indicator 
      if(common.neighbor(x,y,suffStat$skel)){  # 2) x and y have common child
        cond <- TRUE
      }
    }
    return(cond)
  }else{
    return(TRUE)  # 1. if there is skel, we can avoid unnecessary test; else we just do correction.
  }
}

common.neighbor <- function(x,y,skel){
  skel <-as(skel,"matrix")
  sum((skel[,x]==1) & (skel[,y]==1)) > 0  # share the neighbor
}


## ****************** Missing Value PC (MVPC) ******************
# Due to the pMax and sepset we have to change something inside

skeleton2 <- function (suffStat, indepTest, alpha, labels, p,skel_pre, 
                       method = c("stable","original", "stable.fast"), 
                       m.max = Inf, fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE, numCores = 1, verbose = FALSE) 
{
  cl <- match.call()
  if (!missing(p)) 
    stopifnot(is.numeric(p), length(p <- as.integer(p)) == 
                1, p >= 2)
  if (missing(labels)) {
    if (missing(p)) 
      stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  }
  else {
    stopifnot(is.character(labels))
    if (missing(p)) 
      p <- length(labels)
    else if (p != length(labels)) 
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
  }
  seq_p <- seq_len(p)
  method <- match.arg(method)
  if (is.null(fixedGaps)) {
    G <- matrix(TRUE, nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedGaps), c(p, p))) 
    stop("Dimensions of the dataset and fixedGaps do not agree.")
  else if (!identical(fixedGaps, t(fixedGaps))) 
    stop("fixedGaps must be symmetric")
  else G <- !fixedGaps
  diag(G) <- FALSE
  if (any(is.null(fixedEdges))) {
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedEdges), c(p, p))) 
    stop("Dimensions of the dataset and fixedEdges do not agree.")
  else if (!identical(fixedEdges, t(fixedEdges))) 
    stop("fixedEdges must be symmetric")
  stopifnot((is.integer(numCores) || is.numeric(numCores)) && 
              numCores > 0)
  if (numCores > 1 && method != "stable.fast") {
    warning("Argument numCores ignored: parallelization only available for method = 'stable.fast'")
  }
  if (method == "stable.fast") {
    if (identical(indepTest, gaussCItest)) 
      indepTestName <- "gauss"
    else indepTestName <- "rfun"
    options <- list(verbose = as.integer(verbose), m.max = as.integer(ifelse(is.infinite(m.max), 
                                                                             p, m.max)), NAdelete = NAdelete, numCores = numCores)
    res <- .Call("estimateSkeleton", G, suffStat, indepTestName, 
                 indepTest, alpha, fixedEdges, options)
    G <- res$amat
    sepset <- lapply(seq_p, function(i) c(lapply(res$sepset[[i]], 
                                                 function(v) if (identical(v, as.integer(-1))) NULL else v), 
                                          vector("list", p - length(res$sepset[[i]]))))
    pMax <- res$pMax
    n.edgetests <- res$n.edgetests
    ord <- length(n.edgetests) - 1L
  }
  else {
    pval <- NULL
    if(!is.null(skel_pre@sepset)){sepset<-skel_pre@sepset}
    else{sepset <- lapply(seq_p, function(.) vector("list",p))}# a list of lists [p x p]
    
    ## save maximal p value
    if(is.null(skel_pre@sepset)){pMax <- matrix(-Inf, nrow = p, ncol = p)}
    else{pMax <- skel_pre@pMax}
    
    # sepset <- lapply(seq_p, function(.) vector("list", p))
    # pMax <- matrix(-Inf, nrow = p, ncol = p)
    diag(pMax) <- 1
    done <- FALSE
    ord <- 0L
    n.edgetests <- numeric(1)
    while (!done && any(G) && ord <= m.max) {
      n.edgetests[ord1 <- ord + 1L] <- 0
      done <- TRUE
      ind <- which(G, arr.ind = TRUE)
      ind <- ind[order(ind[, 1]), ]
      remEdges <- nrow(ind)
      if (verbose) 
        cat("Order=", ord, "; remaining edges:", remEdges, 
            "\n", sep = "")
      if (method == "stable") {
        G.l <- split(G, gl(p, p))
      }
      for (i in 1:remEdges) {
        if (verbose && (verbose >= 2 || i%%100 == 0)) 
          cat("|i=", i, "|iMax=", remEdges, "\n")
        x <- ind[i, 1]
        y <- ind[i, 2]
        if (G[y, x] && !fixedEdges[y, x]) {
          nbrsBool <- if (method == "stable") 
            G.l[[x]]
          else G[, x]
          nbrsBool[y] <- FALSE
          nbrs <- seq_p[nbrsBool]
          length_nbrs <- length(nbrs)
          if (length_nbrs >= ord) {
            if (length_nbrs > ord) 
              done <- FALSE
            S <- seq_len(ord)
            repeat {
              n.edgetests[ord1] <- n.edgetests[ord1] + 
                1
              pval <- indepTest(x, y, nbrs[S], suffStat)
              if (verbose) 
                cat("x=", x, " y=", y, " S=", nbrs[S], 
                    ": pval =", pval, "\n")
              if (is.na(pval)) 
                pval <- as.numeric(NAdelete)
              if (pMax[x, y] < pval) 
                pMax[x, y] <- pval
              if (pval >= alpha) {
                G[x, y] <- G[y, x] <- FALSE
                sepset[[x]][[y]] <- nbrs[S]
                break
              }
              else {
                nextSet <- getNextSet(length_nbrs, ord, 
                                      S)
                if (nextSet$wasLast) 
                  break
                S <- nextSet$nextSet
              }
            }
          }
        }
      }
      ord <- ord + 1L
    }
    for (i in 1:(p - 1)) {
      for (j in 2:p) pMax[i, j] <- pMax[j, i] <- max(pMax[i,j], pMax[j, i])
    }
  }
  Gobject <- if (sum(G) == 0) {
    new("graphNEL", nodes = labels)
  }
  else {
    colnames(G) <- rownames(G) <- labels
    as(G, "graphNEL")
  }
  new("pcAlgo", graph = Gobject, call = cl, n = integer(0), 
      max.ord = as.integer(ord - 1), n.edgetests = n.edgetests, 
      sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))
}


mvpc<-function(suffStat, indepTest,corrMethod,prt_m, alpha, labels, p,
               fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE, m.max = Inf,
               u2pd = c("relaxed", "rand", "retry"),
               skel.method = c("stable", "original", "stable.fast"),
               conservative = FALSE, maj.rule = FALSE, solve.confl = FALSE, numCores = 1, verbose = FALSE){
  # Test-wise skeleton search result as initialization for detection. 
  # The initialization of skeleton doesnot save that much time for the detection of parents of missingness inidicators. 
  # e.g.:100 000 data points, the time difference is within 1 second.
  # Thus we do not initialize it at the beginning.
  # skel.ini_ <- skeleton(suffStat, indepTest, alpha, p=p)
  # skel.gaps= graph2gaps(skel.ini_)
  
  ## MVPC step1: Detect parents of missingness indicators.
  # prt_m<-get_prt_m_ind(data=suffStat$data, indepTest, alpha, p) # "suffStat$data" is "data_m" which containing missing values.
  # print(prt_m)
  suffStat$prt_m = prt_m
  ## MVPC step2:
  # a) Run PC algorithm with the 1st step skeleton;
  cl <- match.call()
  if (!missing(p))
    stopifnot(is.numeric(p), length(p <- as.integer(p)) ==
                1, p >= 2)
  if (missing(labels)) {
    if (missing(p))
      stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  }
  else {
    stopifnot(is.character(labels))
    if (missing(p)) {
      p <- length(labels)
    }
    else if (p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else message("No need to specify 'p', when 'labels' is given")
  }
  u2pd <- match.arg(u2pd)
  skel.method <- match.arg(skel.method)
  if (u2pd != "relaxed") {
    if (conservative || maj.rule)
      stop("Conservative PC and majority rule PC can only be run with 'u2pd = relaxed'")
    if (solve.confl)
      stop("Versions of PC using lists for the orientation rules (and possibly bi-directed edges)\n can only be run with 'u2pd = relaxed'")
  }
  if (conservative && maj.rule)
    stop("Choose either conservative PC or majority rule PC!")

  skel_pre <- skeleton(suffStat, indepTest, alpha, labels = labels,
                   method = skel.method, fixedGaps = fixedGaps, fixedEdges = fixedEdges,
                   NAdelete = NAdelete, m.max = m.max, numCores = numCores,
                   verbose = verbose)

  suffStat$skel = skel_pre  # For test whether the test need to do the correction
  fixedGaps <- graph2gaps(skel_pre)

  # b) Correction of the extra edges
  
  skel <- skeleton2(suffStat, corrMethod, alpha, skel_pre=skel_pre, labels = labels,
                   method = skel.method, fixedGaps = fixedGaps, fixedEdges = fixedEdges,
                   NAdelete = NAdelete, m.max = m.max, numCores = numCores,
                   verbose = verbose)
  
  # c) Orient the edges
  skel@call <- cl
  if (!conservative && !maj.rule) {
    switch(u2pd, rand = udag2pdag(skel), retry = udag2pdagSpecial(skel)$pcObj,
           relaxed = udag2pdagRelaxed(skel, verbose = verbose,
                                      solve.confl = solve.confl))
  }
  else {
    pc. <- pc.cons.intern(skel, suffStat, indepTest, alpha,
                          version.unf = c(2, 1), maj.rule = maj.rule, verbose = verbose)
    udag2pdagRelaxed(pc.$sk, verbose = verbose, unfVect = pc.$unfTripl,
                     solve.confl = solve.confl)
  }
}

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

graph2gaps <- function(graphnet){
  # graphnet: a graphNet class which is defined in pcalg
  #
  # ******************
  # Return: a gaps matrix which represent the gaps in a DAG or a skeleton of a graph
  graphnet.matrix <-as(graphnet@graph,"matrix")
  gaps <- graphnet.matrix ==0
  return(gaps)
}

get_prt_m_ind <- function(data, indepTest, alpha, p, fixedGaps = NULL){
  ## Detect causes of the missingness indicators
  ## data:
  ## mode: the best is Anova test to test the conditional independence (the test for continuous variables also works)
  ## ---------------------------
  ## Return: the index of parents of missingness indicators
  R<-get_m_ind(data)
  m=c()
  prt=list()
  suffStat = list(data=data)
  
  count=1
  for(R_ind in R){
    prt_R_ind<- get_prt_R_ind(suffStat, indepTest, alpha, p, R_ind, fixedGaps=fixedGaps)
    if(length(prt_R_ind)!=0){
      m[count] <- R_ind
      prt[[count]]<- prt_R_ind
      count = count + 1 
    }
  }
  prt_m<-data.frame(m=m)
  prt_m[['prt']]<-prt
  return(prt_m)
}

get_prt_R_ind <- function(suffStat, indepTest, alpha, p, R_ind,
                          labels,  method = c("stable", "original", "stable.fast"), 
                          m.max = Inf, fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE, numCores = 1, verbose = FALSE){
  # the function can be regared as "<- function(suffStat, indepTest, alpha, p, R_ind, skel.ini)"
  # This is the skeleton search algorithm adapted from skeleton(...) in pcalg
  # 1. Change data of the R_ind column into binary missingness indicator, and the corresponding graph
  # 2. Only test the edges between R_ind and all the other substantial variables 
  # 3. Test with the Anovatest (because all the tests consist of a binary variable and a continuous variable give a set of continuous variables)
  
  cl <- match.call()
  if (!missing(p)) 
    stopifnot(is.numeric(p), length(p <- as.integer(p)) == 
                1, p >= 2)
  if (missing(labels)) {
    if (missing(p)) 
      stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  }
  else {
    stopifnot(is.character(labels))
    if (missing(p)) 
      p <- length(labels)
    else if (p != length(labels)) 
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
  }
  seq_p <- seq_len(p)
  method <- match.arg(method)
  if (is.null(fixedGaps)) {
    G <- matrix(TRUE, nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedGaps), c(p, p))) 
    stop("Dimensions of the dataset and fixedGaps do not agree.")
  else if (!identical(fixedGaps, t(fixedGaps))) 
    stop("fixedGaps must be symmetric")
  else G <- !fixedGaps
  
  # ****** Begin adapted 1 ******
  G[,R_ind]<-TRUE
  G[R_ind,]<-TRUE
  suffStat$data[,R_ind]<-as.integer(is.na(suffStat$data[,R_ind]))
  # ****** End adapted 1 ******
  diag(G) <- FALSE
  
  if (any(is.null(fixedEdges))) {
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedEdges), c(p, p))) 
    stop("Dimensions of the dataset and fixedEdges do not agree.")
  else if (!identical(fixedEdges, t(fixedEdges))) 
    stop("fixedEdges must be symmetric")
  stopifnot((is.integer(numCores) || is.numeric(numCores)) && 
              numCores > 0)
  if (numCores > 1 && method != "stable.fast") {
    warning("Argument numCores ignored: parallelization only available for method = 'stable.fast'")
  }
  if (method == "stable.fast") {
    if (identical(indepTest, gaussCItest)) 
      indepTestName <- "gauss"
    else indepTestName <- "rfun"
    options <- list(verbose = as.integer(verbose), m.max = as.integer(ifelse(is.infinite(m.max), 
                                                                             p, m.max)), NAdelete = NAdelete, numCores = numCores)
    res <- .Call("estimateSkeleton", G, suffStat, indepTestName, 
                 indepTest, alpha, fixedEdges, options)
    G <- res$amat
    sepset <- lapply(seq_p, function(i) c(lapply(res$sepset[[i]], 
                                                 function(v) if (identical(v, as.integer(-1))) NULL else v), 
                                          vector("list", p - length(res$sepset[[i]]))))
    pMax <- res$pMax
    n.edgetests <- res$n.edgetests
    ord <- length(n.edgetests) - 1L
  }
  else {
    pval <- NULL
    sepset <- lapply(seq_p, function(.) vector("list", p))
    pMax <- matrix(-Inf, nrow = p, ncol = p)
    diag(pMax) <- 1
    done <- FALSE
    ord <- 0L
    n.edgetests <- numeric(1)
    while (!done && any(G) && ord <= m.max) {
      n.edgetests[ord1 <- ord + 1L] <- 0
      done <- TRUE
      ind <- which(G, arr.ind = TRUE)
      ind <- ind[order(ind[, 1]), ]
      remEdges <- nrow(ind)
      if (verbose) 
        cat("Order=", ord, "; remaining edges:", remEdges, 
            "\n", sep = "")
      if (method == "stable") {
        G.l <- split(G, gl(p, p))
      }
      for (i in 1:remEdges) {
        if (verbose && (verbose >= 2 || i%%100 == 0)) 
          cat("|i=", i, "|iMax=", remEdges, "\n")
        x <- ind[i, 1]
        y <- ind[i, 2]
        
        # ****** Begin adapted 2 ******
        if(x!=R_ind){next}  # we only focus on the missingness variable
        # ****** End adapted 2 ******
        
        if (G[y, x] && !fixedEdges[y, x]) {
          nbrsBool <- if (method == "stable") 
            G.l[[x]]
          else G[, x]
          nbrsBool[y] <- FALSE
          nbrs <- seq_p[nbrsBool]
          length_nbrs <- length(nbrs)
          if (length_nbrs >= ord) {
            if (length_nbrs > ord) 
              done <- FALSE
            S <- seq_len(ord)
            repeat {
              n.edgetests[ord1] <- n.edgetests[ord1] + 
                1
              pval <- indepTest(x, y, nbrs[S], suffStat)
              if (verbose) 
                cat("x=", x, " y=", y, " S=", nbrs[S], 
                    ": pval =", pval, "\n")
              if (is.na(pval)) 
                pval <- as.numeric(NAdelete)
              if (pMax[x, y] < pval) 
                pMax[x, y] <- pval
              if (pval >= alpha) {
                G[x, y] <- G[y, x] <- FALSE
                sepset[[x]][[y]] <- nbrs[S]
                break
              }
              else {
                nextSet <- getNextSet(length_nbrs, ord, 
                                      S)
                if (nextSet$wasLast) 
                  break
                S <- nextSet$nextSet
              }
            }
          }
        }
      }
      ord <- ord + 1L
    }
    for (i in 1:(p - 1)) {
      for (j in 2:p) pMax[i, j] <- pMax[j, i] <- max(pMax[i, 
                                                          j], pMax[j, i])
    }
  }
  
  prts<-G[R_ind,]
  x<-c(1:length(suffStat$data[1,]))
  return(x[prts])
}


#****************** Evaluation ******************
comp_com_td_mvpc<- function(exp="mar_cor",
                            n=1000, 
                            times=10,
                            p=20,
                            data_path="/Users/ruibo/Desktop/mvpc/mvpc-xyz/causality/R-proj/data"
                            ){
  shd_comp= rep(0,times)
  shd_mvpc = rep(0,times)
  shd_test = rep(0,times)
  
  for(i in 1:times){
    print(i)
    # i=10
    load(paste(data_path,"/syn/",exp,n,"_",i,".rda",sep=""))
    
    data_mcar = data_all$data_mcar
    data_c = data_all$data_c
    suffStat = list(data_m = data_all$suffStat$data)
    myCPDAG=data_all$cpdag
    myDAG=data_all$dag
    
    res_mvpc<-mvpc(suffStat, gaussCItest_td,PermCCItest, alpha=0.01, p=p)
    res_pc<-pc(suffStat, gaussCItest_td, alpha=0.01, p=p)
    
    suffStat  = list(C = cor(data_c),n=length(data_c[,1]))
    res_ref<-pc(suffStat, gaussCItest, alpha=0.01, p=p)
    
    shd_comp[i] = shd(res_ref, myCPDAG)
    shd_mvpc[i] = shd(res_mvpc, myCPDAG)
    shd_test[i] = shd(res_pc, myCPDAG)  
  }
  return(list(shd_comp=shd_comp,shd_mvpc=shd_mvpc,shd_test=shd_test))
}

eva.detection<-function(prt1,prt2){
  count = 0
  for(i in 1:length(prt1)){
    if(length(prt1[i])!=length(prt2[[i]])){count = count + 1 }
    else if(prt1[i]!=prt2[[i]]){count = count + 1 }
  }
  return(count)
}

test_adj<-function(myCPDAG,res){
  ## Ajacency
  
  true_skel<-get_edge_pairs(myCPDAG)
  our_skel<-get_edge_pairs(res)
  
  tran_true_skel<-get_trans_edge_pairs(myCPDAG)
  tran_our_skel<-get_trans_edge_pairs(res)
  
  true_pair <- complex(real = true_skel[,1], imaginary = true_skel[,2])
  ours_pair <- complex(real = our_skel[,1], imaginary = our_skel[,2])
  
  tran_true_pair <- complex(real = tran_true_skel[,1], imaginary = tran_true_skel[,2])
  tran_ours_pair <- complex(real = tran_our_skel[,1], imaginary = tran_our_skel[,2])
  
  # all pairs:
  all_pair<-union(true_pair,tran_true_pair)
  all_pair_us <-union(ours_pair,tran_ours_pair)
  
  num_our<-length(all_pair_us)
  num_true<-length(all_pair)
  num_us_correct<-length(intersect(all_pair_us,all_pair))
  recall<-num_us_correct/num_true
  precision<-num_us_correct/num_our
  return(list(recall=recall, precision=precision))
}

get_edge_pairs<-function(G){
  g<-as(G,'matrix')
  g[lower.tri(g)]<-0
  g<-g==1
  which(g, arr.ind = TRUE)
}

get_trans_edge_pairs<-function(G){
  g<-t(as(G,'matrix'))
  g[lower.tri(g)]<-0
  g<-g==1
  which(g, arr.ind = TRUE)
}

compute_rp<-function(rp_){
  recall = c()
  precision = c()
  count = 1
  for(i in 1:length(rp_)){
    if(i == 6 | i ==12 | i==18 | i==24){next}
    recall[count] = rp_[[i]]$recall
    precision[count] = rp_[[i]]$precision
    count = count + 1
  }
  return(c(mean(recall), sd(recall), mean(precision),sd(precision)))
}

