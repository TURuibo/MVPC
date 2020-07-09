library(pcalg)
library(mvtnorm)
library(weights)
library(ks)
library(e1071)  # For generating all combinations of a binary vector
library(weights)  # For kernel density estimate
library(ks)
library(DescTools)
library(mipfp)

## ****************** The famework MVPC: Missing Value PC (MVPC) ******************
mvpc<-function(suffStat, indepTest,corrMethod,alpha,p,
               prt_m = NULL,
               labels=NULL,
               fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE, m.max = Inf,
               u2pd = c("relaxed", "rand", "retry"),
               skel.method = c("stable", "original", "stable.fast"),
               conservative = FALSE, maj.rule = FALSE, solve.confl = FALSE, numCores = 1, verbose = FALSE){
  # Input:
  #   suffStat: 1) a list of sufficient statitics for the conditional independence tests such as gaussCItest
  #             2) a dataset, named data for the conditional independence tests, 
  #                 gaussCItest.td, gaussCItest.drw, gaussCItest.permc, binCItest.permc, binCItest.drw
  #   indepTest: conditional independence test:
  #               gaussCItest, gaussCItest.td
  #   corrMethod: correction methods for continuous variable: gaussCItest.drw, gaussCItest.permc,
  #               correction methods for binary variable: binCItest.permc, binCItest.drw
  #   alpha: p-value for conditional independence tests 
  #   p: the number of variables
  #   prt_m: a data.frame, which saves the missingness indicators in prt_m$m and their parent in prt_m$prt.
  #     prt_m$m: a list of the missingness indicators.
  #     prt_m$prt: a collection (list) of lists which are the parents of corresponding missingness indicators
  # Other inputs and Return follow the function pc in the package "pcalg"
  
  ## MVPC step1: Detect parents of missingness indicators.
  prt_m<-detection_prt_m(data=suffStat$data, indepTest, alpha, p) # "suffStat$data" is "data_m" which containing missing values.
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

## Detection of the missingness indicator parents
detection_prt_m <- function(data, indepTest, alpha, p, fixedGaps = NULL){
  ## Return: the missingness indicator and their parents as a dataframe prt_m for mvpc
  #   prt_m: a data.frame, which saves the missingness indicators in prt_m$m and their parent in prt_m$prt.
  #     prt_m$m: a list of the missingness indicators.
  #     prt_m$prt: a collection (list) of lists which are the parents of corresponding missingness indicators
  
  return(get_prt_m_ind(data, indepTest, alpha, p, fixedGaps))
}

## Causal skeleton search with correction methods
skeleton2 <- function (suffStat, indepTest, alpha, labels, p,skel_pre, 
                       method = c("stable","original", "stable.fast"), 
                       m.max = Inf, fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE, numCores = 1, verbose = FALSE){
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
  ## Return: the missingness indicator and their parents as a dataframe prt_m for mvpc
  #   prt_m: a data.frame, which saves the missingness indicators in prt_m$m and their parent in prt_m$prt.
  #     prt_m$m: a list of the missingness indicators.
  #     prt_m$prt: a collection (list) of lists which are the parents of corresponding missingness indicators
  
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


#************ set up the environment of R.Matlab ************
initialize.rmatlab<-function(){
  library(R.matlab)
  print("Setting up matlab")
  options(matlab="/Applications/MATLAB_R2019b.app/bin/matlab")
  Matlab$startServer()
  matlab <- Matlab()
  isOpen <- open(matlab)
  path = "/Users/ruibo/Desktop/mvpc/mvpc_v2/MVPC/src/KCI-test"
  addpath = paste0("addpath('",path,"')", sep = "")
  evaluate(matlab, addpath)
  return(matlab)
}
# Example: 
# matlab <- initialize.rmatlab()
# close(matlab)
# #************ END: set up the environment of R.Matlab ************  