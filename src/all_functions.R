library(pcalg)
library(Rfast)

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
  
  myDAG <- randomDAG(p,2/(p-1))
  ## Make sure the graph contains collider
  # while(!check_collider(myDAG)){
  #   myDAG <- randomDAG(p,2/(p-1))
  # }
  set.seed(seed)
  data_complete <- rmvnorm(n, mu=rep(0,p), sigma=trueCov(myDAG))
  
  # Make sure the collider is the cause of R-variable of its parent variable
  data_del <-del(data_complete,myDAG,mode)
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

del <- function(data, myDAG,mode='mar',num_m=6){
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
  
  size <- dim(data)
  num_var <- size[2]
  num_samples <- size[1]
  # 1. Choose the missingness indicators/ choose which variables contain missing values 
  m_ind = select_m_ind(num_var,num_m)
  
  # 2. Choose the cause of missingess indicators
  if(mode=="mar"){
    prt_m = select_prt_mar_ind(m_ind,num_var)
  }else{
    prt_m = select_prt_mnar_ind(m_ind,num_var)
  }
  
  # 3. Generate missing values according the the values of missingness indicators
  data_mar_ref = generate_missing_values(data,m_ind,prt_m)
  data_mar = data_mar_ref$data_mar
  data_ref = data_mar_ref$data_ref
  
  return(list(data_mar=data_mar,
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

generate_missing_values<-function(data,m_ind,parent_m_ind){
  # Give the parents of missingness indicators, and the missingness indcators
  # The missing values are generated whn the parent values are in the bottom XX percentage.
  
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

# Example of generating synthetic data
# gen_data(20,10,"mnar")

#******************Functions of extracting necessary information for MVPC******************
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

get_prt_m_ind<-function(R,data,alpha=0.01,mode=Anovatest){
  ## Detect causes of the missingness indicators
  ## R: 
  ## data:
  ## mode: Anova test to test the conditional independence
  ## ---------------------------
  ## Return: the index of parents of missingness indicators

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


skeleton3 <- function(suffStat,alpha=0.05, R_ind,labels,
                      m_method=deletion,
                      indepTest=Anovatest,
                      method = c("stable", "original", "stable.fast"), m.max = Inf,
                      fixedGaps = NULL, fixedEdges = NULL,
                      NAdelete = TRUE, numCores = 1, verbose = FALSE,debug=FALSE)
{
  ## Assumption : 1. R varaible can only be the child rather than parent variable
  ##            : 2. R only has one cause at most.
  ## Purpose: Find the variables that are connected to a certain R varaible
  ## estimate skeleton of DAG given data
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - suffStat: List containing all necessary elements for the conditional
  ##             independence decisions in the function "indepTest".
  ## - indepTest: predefined function for testing conditional independence
  ## - alpha: Significance level of individual partial correlation tests
  ## ----------------------------------------------------------------------
  ## Return: index of the support variables (the cause of the missingness indicator)
  
  r <- is.na(suffStat$data[,R_ind])
  suffStat$data[,R_ind]<-rep(0,length(suffStat$data[,R_ind]))
  suffStat$data[r,R_ind]<-1
  suffStat$C = cor(suffStat$data)
  m_var<-get_m_ind(suffStat$data)
  p<-length(suffStat$data[1,])
  
  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p))
      p <- length(labels)
    else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    ## Don't want message, in case this is called e.g. from fciPlus():
    ## else
    ##   message("No need to specify 'p', when 'labels' is given")
  }
  seq_p <- seq_len(p)
  method <- match.arg(method)
  ## C++ version still has problems under Windows; will have to check why
  ##  if (method == "stable.fast" && .Platform$OS.type == "windows") {
  ##    method <- "stable"
  ##    warning("Method 'stable.fast' is not available under Windows; using 'stable' instead.")
  ##  }
  
  ## G := !fixedGaps, i.e. G[i,j] is true  iff  i--j  will be investigated
  if (is.null(fixedGaps)) {
    ##
    G <- matrix(TRUE, nrow = p, ncol = p)
    ##
  }
  else if (!identical(dim(fixedGaps), c(p, p)))
    stop("Dimensions of the dataset and fixedGaps do not agree.")
  else if (!identical(fixedGaps, t(fixedGaps)) )
    stop("fixedGaps must be symmetric")
  else
    G <- !fixedGaps
  ## Make sure G is a symetric matrix
  G[,R_ind]<-TRUE
  G[R_ind,]<-TRUE
  diag(G) <- FALSE
  if (any(is.null(fixedEdges))) { ## MM: could be sparse
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedEdges), c(p, p)))
    stop("Dimensions of the dataset and fixedEdges do not agree.")
  else if (!identical(fixedEdges, t(fixedEdges)) )
    stop("fixedEdges must be symmetric")
  
  ## Check number of cores
  stopifnot((is.integer(numCores) || is.numeric(numCores)) && numCores > 0)
  if (numCores > 1 && method != "stable.fast") {
    warning("Argument numCores ignored: parallelization only available for method = 'stable.fast'")
  }
  if (method == "stable.fast") {
    ## Do calculation in C++...
    next
  }
  else {
    ## Original R version
    
    pval <- NULL
    sepset <- lapply(seq_p, function(.) vector("list",p))# a list of lists [p x p]
    ## save maximal p value
    pMax <- matrix(-Inf, nrow = p, ncol = p)
    diag(pMax) <- 1
    done <- FALSE
    ord <- 0L
    n.edgetests <- numeric(1)# final length = max { ord}
    while (!done && any(G) && ord <= m.max) {
      n.edgetests[ord1 <- ord+1L] <- 0
      done <- TRUE
      ind <- which(G, arr.ind = TRUE)
      ## For comparison with C++ sort according to first row
      ind <- ind[order(ind[, 1]), ]
      remEdges <- nrow(ind)
      if (verbose)
        cat("Order=", ord, "; remaining edges:", remEdges,"\n",sep = "")
      if(method == "stable") {
        ## Order-independent version: Compute the adjacency sets for any vertex
        ## Then don't update when edges are deleted
        G.l <- split(G, gl(p,p))
      }
      for (i in 1:remEdges) {
        if(verbose && (verbose >= 2 || i%%100 == 0)) cat("|i=", i, "|iMax=", remEdges, "\n")
        x <- ind[i, 1]
        y <- ind[i, 2]
        if(x!=R_ind){next}  # we only focus on the missingness variable
        if (G[y, x] && !fixedEdges[y, x]) {
          nbrsBool <- if(method == "stable") G.l[[x]] else G[,x]
          nbrsBool[y] <- FALSE
          nbrs <- seq_p[nbrsBool]
          length_nbrs <- length(nbrs)
          if (length_nbrs >= ord) {
            if (length_nbrs > ord)
              done <- FALSE
            S <- seq_len(ord)
            repeat { ## condition w.r.to all  nbrs[S] of size 'ord'
              n.edgetests[ord1] <- n.edgetests[ord1] + 1
              if (identical(indepTest, gaussCItest)){
                pval <- indepTest(x, y, nbrs[S], suffStat) 
              }
              else{
                pval <- indepTest(x, y, nbrs[S], suffStat,m_method)
              }
              if (verbose)
                cat("x=", x, " y=", y, " S=", nbrs[S], ": pval =", pval, "\n")
              if(is.na(pval))
                pval <- as.numeric(NAdelete) ## = if(NAdelete) 1 else 0
              if (pMax[x, y] < pval)
                pMax[x, y] <- pval
              if(identical(indepTest, BICtest)||identical(indepTest, DICtest)){
                if(pval >= 0) { # independent
                  G[x, y] <- G[y, x] <- FALSE
                  sepset[[x]][[y]] <- nbrs[S]
                  
                  break
                }
                else {
                  nextSet <- getNextSet(length_nbrs, ord, S)
                  if (nextSet$wasLast)
                    break
                  S <- nextSet$nextSet
                }
              }
              else{
                if(pval >= alpha) { # independent
                  G[x, y] <- G[y, x] <- FALSE
                  sepset[[x]][[y]] <- nbrs[S]
                  break
                }
                else {
                  nextSet <- getNextSet(length_nbrs, ord, S)
                  if (nextSet$wasLast)
                    break
                  S <- nextSet$nextSet
                }
              }
              
            } ## repeat
          }
        }
      }# for( i )
      ord <- ord + 1L
    } ## while()
  }
  sup<-G[R_ind,]
  x<-c(1:length(suffStat$data[1,]))
  x[sup]
}## end{ skeleton3}


Anovatest <- function(x,y,S,suffStat,method=deletion) {
  ## Conditional independence test (ANOVA) between continuous variables with deletion methods
  ## test P(x,y|S)
  ## suffStat: the class contains the dataset and other necessary variables
  ##    suffStat$C: correlation matrix
  ##    suffStat$n: sample size
  ##    suffStat$data: the dataset
  ##--------------
  ## Return: the p-value of the test 
  
  t_var = c(x,y,S)
  
  if(check_m_var(x,y,S,suffStat)){data = method(t_var,suffStat)}
  else{data<-suffStat$data[,c(x,y,S)]}
  
  if(length(t_var)>2){
    xnam <- paste0("data[,", 3:length(t_var),"]")
    fmla <- as.formula(paste("data[,2] ~ ", paste(xnam, collapse= "+"),"+data[,1]"))
    data[,1]<-factor(data[,1])
    if(length(levels(data[,1]))<=1){
      0
    }else{
      result <- aov(fmla,data = data)
      summary(result)[[1]][length(t_var)-1,5]
    }
    
  }
  else{
    fmla <- as.formula(paste("data[,2] ~ data[,1]"))
    data[,1]<-factor(data[,1])
    
    if(length(levels(data[,1]))<=1){
      1
    }
    else if(sum(data[,1]==levels(data[,1])[1])>0.95*length(suffStat$data[,1])
            |sum(data[,1]==levels(data[,1])[2])>0.95*length(suffStat$data[,1])
            |sum(data[,1]==levels(data[,1])[1])<0.05*length(suffStat$data[,1])
            |sum(data[,1]==levels(data[,1])[2])<0.05*length(suffStat$data[,1])
    ){0}
    else{
      if(length(unique(data[,2]))==2){
        v<-unique(data[,2])
        data<-as(data,'matrix')
        data<-apply(data, 2, as.numeric)
        v<-unique(data[,2])
        if(sum(data[,2]==v[1])>0.9*length(data[,1])|sum(data[,2]==v[2])>0.9*length(data[,1]) |sum(data[,2]==v[1])<0.1*length(data[,1])|sum(data[,2]==v[2])<0.1*length(data[,1])){
          0
        }else{
          binCItest(1,2,c(), list(dm = data, adaptDF = FALSE))
        }
      }else{
        v<-unique(data[,2])
        
        if(sum(data[,2]==v[1])>0.9*length(data[,1])|sum(data[,2]==v[2])>0.9*length(data[,1]) |sum(data[,2]==v[1])<0.1*length(data[,1])|sum(data[,2]==v[2])<0.1*length(data[,1])){
          0
        }else{
          ttest=t.test(fmla,data = data)
          ttest[['p.value']]
        }
      }
      
    }
  }
}