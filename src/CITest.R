#****************** (Conditional) Independence Test ****************** 

gaussCItest.td.ref <- function(x, y, S, suffStat) {
  ## Conditional independence test between continuous variables with deletion methods
  ## test P(x,y|S)
  ## suffStat: the class contains the dataset and other necessary variables
  ##    suffStat$data: the dataset
  ##--------------
  ## Return: the p-value of the test 
  
  data = test_wise_deletion(c(x,y,S),suffStat$data)
  sample_size <<- c(sample_size,length(data[,1]))
  
  suffStat$data = data
  suffStat$C = cor(data)
  suffStat$n = length(data[,1])
  gaussCItest(x, y, S, suffStat)
}

gaussCItest.td <- function(x, y, S, suffStat) {
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

gaussCItest.permc <- function(x, y, S, suffStat){
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
    pval = gaussCItest(1, 2, 3:length(ind_test), suffStat_perm)
    pval
  }
  else{
    pval = gaussCItest(1, 2,c(), suffStat_perm)
    pval
  }
}

binCItest.td.ref<- function(x, y, S, suffStat){
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

binCItest.td <- function(x, y, S, suffStat){
  test_ind = c(x,y,S)
  data = test_wise_deletion(test_ind,suffStat$data)
  if(length(S)>0){
    binCItest(1, 2, 3:length(test_ind), list(dm = data[,test_ind], adaptDF = FALSE))  
  }
  else{
    binCItest(1, 2, c(), list(dm = data[,test_ind], adaptDF = FALSE))
  }
  
}

gaussCItest.drw <- function(x, y, S, suffStat) {
  ## Conditional independence test between continuous variables with deletion methods
  ## test P(x,y|S)
  ## suffStat: the class contains the dataset and other necessary variables
  ##    suffStat$data: the dataset
  ##--------------
  ## Return: the p-value of the test 
  # Note that "tw_data" and "weights" should have the same order/index
  if(!cond.PermC(x, y, S, suffStat)){return(gaussCItest_td(x,y,S,suffStat))}
  
  ind_W <- unique(get_prt_m_xys(c(x,y,S), suffStat))  # Get parents the {xyS} missingness indicators
  if(length(ind_W)==0){return(gaussCItest_td(x,y,S,suffStat))}
  
  pa_W <- unique(get_prt_m_xys(ind_W, suffStat))
  candi_W <- setdiff(pa_W, ind_W)
  
  while(length(candi_W) > 0  ){
    ind_W <- c(ind_W, candi_W) # Get parents the W missingness indicators
    pa_W <- unique(get_prt_m_xys(ind_W, suffStat))
    candi_W <- setdiff(pa_W, ind_W)
  } 
  ind_W <- unique(ind_W)
  
  corr_ind = c(x,y,S,ind_W)
  tw_data = test_wise_deletion(corr_ind, suffStat$data) # Fixed bug: The same number of variables as the orginal dataset
  weights = compute.weights.continuous(corr_ind, suffStat)
  
  suffStat$C = wtd.cors(tw_data, tw_data, weights)
  suffStat$n = length(weights)
  gaussCItest(x, y, S,suffStat)
}

binCItest.permc<- function(x, y, S, suffStat){  
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
      if(length(S) > 0){
        pval = binCItest(1,2,c(3:length(ind_test)), list(dm = data.vir, adaptDF = TRUE))
        pval
      }
      else{
        pval = binCItest(1,2,c(), list(dm = data.vir , adaptDF = TRUE))
        pval
      }
    },
    error=function(cond) {
      message(cond)
      # Choose a return value in case of error
      return(binCItest_td(x,y,S,suffStat))
    })
  
}  

binCItest.drw <- function(x, y, S, suffStat){
  if(!cond.PermC(x, y, S, suffStat)){return(binCItest_td(x,y,S,suffStat))}
  weights_ = compute_weights(x, y, S, suffStat)
  test_ind = c(x,y,S)
  del_res = test_wise_deletion_w(test_ind,suffStat$data, weights_)
  
  data = del_res$data
  weights = del_res$weights
  if(length(S)>0){
    pval = binCItest_w(1, 2, 3:length(test_ind), weights,list(dm = data[,test_ind], adaptDF = FALSE))
    pval
    
  }
  else{
    pval = binCItest_w(1, 2, c(), weights, list(dm = data[,test_ind], adaptDF = FALSE))
    pval
  }
}


#****************** Util for the CI tests ****************** 

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
  for(ri in unique(c(x,y,S,ind_W))){
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
  
  pa_W <- unique(get_prt_m_xys(ind_W, suffStat))
  candi_W <- setdiff(pa_W, ind_W)
  
  while(length(candi_W) > 0  ){
    ind_W <- c(ind_W, candi_W) # Get parents the W missingness indicators
    pa_W <- unique(get_prt_m_xys(ind_W, suffStat))
    candi_W <- setdiff(pa_W, ind_W)
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
    weights[is.na(ind)] = NA
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

compute.weights.continuous<-function(corr_ind, suffStat, kernel.method = kde.weights){
  ind.twdel = indx_test_wise_deletion(corr_ind,suffStat$data)
  weights = rep(1,length(ind.twdel))  # length of test-wise deleted data
  ind_r_xys = get_ind_r_xys(corr_ind, suffStat)
  
  for(ind_ri in ind_r_xys){
    ind_pa = get_prt_i(ind_ri,suffStat)
    # ind_pa =get_prt_m_xys(c(ind_ri), suffStat) 
    # Return the single parent of a missingness indicator (an element, i.e., list[[i]], not a list, i.e., list[i])
    pa = suffStat$data[,ind_pa]
    beta = kernel.method(pa[ind.twdel], pa[!is.na(pa)])
    weights = weights * beta
  }
  weights
}

kde.weights <- function(x_del,x_full){
  f_w= approxfun(density(x_full))
  f_wr= approxfun(density(x_del))
  beta=f_w(x_del)/f_wr(x_del)
  beta
}

kdrw.weights <- function (x, x_target){
  # kernel-based density ratio estimate
  setVariable(matlab, X = x)
  setVariable(matlab, X_target=x_target)
  if(length(x)<800){
    coeff = 0.7
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
  beta = beta
}

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
