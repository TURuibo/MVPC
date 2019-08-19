library(pcalg)
library(mvtnorm)
library(weights)
library(ks)


#******************Functions for Synthetic Data Generation******************

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

select_m_ind <- function(num_var,num_m){
  return(sample(x=1:num_var,size=num_m,replace=FALSE))
}

select_prt_mar_ind <- function(m_ind,num_var){
  return(sample(x=setdiff(1:num_var,m_ind),size=length(m_ind),replace=TRUE))
}

select_prt_mnar_ind <- function(m_ind,num_var){
  # No self-masking
  # only the variable contains missing values can be the parents of missingness indicators
  prts = c()
  for(i in m_ind){
    newelem<-sample(setdiff(m_ind,i),size=1)
    prts <- c(prts, newelem) 
  }
  return(prts)
}

generate_missing_values <- function(p,n,myDAG,m_ind,parent_m_ind){
  # Give the parents of missingness indicators, and the missingness indcators
  # The missing values are generated whn the parent values are in the bottom XX percentage.
  data <- rmvnorm(n, mean=rep(0,p), sigma=trueCov(myDAG))
  data_m=data
  data_mcar=data
  for(i in c(1:length(m_ind))){
    # Choose lower "bottom_p" percentage of the values
    bottom_p <- runif(1, min = 0.1, max = 0.5)
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

# Example of generating synthetic data
# gen_data(20,10,"mnar")


#****************** (Conditional) Independence Test ****************** 


binCItest_td <- function(x, y, S, suffStat){
  data = test_wise_deletion(c(x,y,S),suffStat$dm)
  suffStat$dm = data
  binCItest(x, y, S, suffStat)
}

DRWbinCItest <- function(x, y, S, suffStat){
  
  
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
  # z <- zStat(x,y,S, 
  #            C = cor(data), n = length(data[,1]))
  # 2*pnorm(abs(z), lower.tail = FALSE)
}

PermCCItest <- function(x, y, S, suffStat){
  ## The Z <- XY; Rz <- XY is not included in the test 
  ## Step 1: Learning generaive model for {X, Y, S} to impute X, Y, and S
  W = get_prt_m_xys(c(x,y,S), suffStat)  # Get parents the {xyS} missingness indicators: prt_m
  ind_permc <- c(x, y, S, W)
  ind_test <- c(x, y, S)
  data <- test_wise_deletion(ind_permc, suffStat$data)
  data <- data[,ind_permc]
  
  xnam <- paste0("data[,", (length(ind_test)+1):length(ind_permc),"]")
  lr <- list()
  res <- list()
  for( i in ind_test){
    fmla <- as.formula(paste("data[,", i,"] ~ 0 + ", paste(xnam, collapse= "+"),""))
    lr[[i]] <- lm(formula = fmla, data = data)  # the ith block of list is ...
    res[[i]] <- residuals(lr[[i]])  # residuals
  }
  
  ## Step 2: Shuffle the source "W" -- parents of the missingness indicators
  #         a) Remove the missing entries of dat[, W]
  #         b) Row-based/ sample based Shuffle the index of W data points
  #         c) the same number of test-wise deletion CI Test 
  data_W_p = perm(W, suffStat$data)
  data[,(length(ind_test)+1):length(ind_permc)] = data_W_p[1:length(data[,1]),]
  
  ## Step 3: Generate the virtual data follows the full data distribution P(X, Y, S)
  vir <- list()
  for(i in ind_test){
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
  if(! iscorr(x, y, S, suffStat)){
    gaussCItest_td(x, y, S, suffStat)
  }
  else{
    # Logistic regression for each missingness indicator (rx, ry, rsi): ri ~ Pa(ri).
    # Note that data for training the logistic regression are different with the data used for computing weights.
    ind_r_xys <- get_ind_r_xys(x, y, S, suffStat)
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

get_ind_r_xys <- function(x, y, S, suffStat){
  # return index of missingness indicators of x, y, S
  ind_r_xys <- c()
  for(i in c(x,y,S)){
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
#****************** Missing Value PC (MVPC) ******************

mvpc<-function(suffStat, indepTest, alpha, labels, p,
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
  prt_m<-get_prt_m_ind(data=suffStat$data, indepTest, alpha, p) # "suffStat$data" is "data_m" which containing missing values.
  
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

  skel <- skeleton(suffStat, indepTest, alpha, labels = labels,
                   method = skel.method, fixedGaps = fixedGaps, fixedEdges = fixedEdges,
                   NAdelete = NAdelete, m.max = m.max, numCores = numCores,
                   verbose = verbose)

  suffStat$skel = skel  # For test whether the test need to do the correction
  fixedGaps <- graph2gaps(skel)

  # b) Correction of the extra edges
  skel <- skeleton(suffStat, indepTest, alpha, labels = labels,
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
  R<-get_m_ind(data_m)
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

eva.detection<-function(prt1,prt2){
  count = 0
  for(i in 1:length(prt1)){
    if(length(prt1[i])!=length(prt2[[i]])){count = count + 1 }
    else if(prt1[i]!=prt2[[i]]){count = count + 1 }
  }
  return(count)
}
