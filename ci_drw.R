# binCItest_drw <- function (x, y, S, dm, c0,c1,adaptDF = FALSE, n.min = 10 * df, verbose = FALSE) 
# {
#   stopifnot((n <- nrow(dm)) >= 1)
#   if (!all(as.logical(dm) == dm)) 
#     stop("'dm' must be binary, i.e. with values in {0,1}")
#   if (verbose) 
#     cat("Edge ", x, "--", y, " with subset S =", S, "\n")
#   lenS <- length(S)
#   df <- 2^lenS
#   if (n < n.min) {
#     warning(gettextf("n=%d is too small (n < n.min = %d ) for G^2 test (=> treated as independence)", 
#                      n, n.min), domain = NA)
#     return(1)
#   }
#   d.x1 <- dm[, x] + 1L
#   d.y1 <- dm[, y] + 1L
#   if (lenS <= 5) {
#     n12 <- 1:2
#     switch(lenS + 1L, {
#       nijk = cpd_1(d.x1[w==0], d.y1[w==0]) * c0 + cpd_1(d.x1[w==1], d.y1[w==1]) * c1
#       # nijk <- array(0L, c(2, 2))
#       # for (i in n12) {
#       #   d.x.i <- d.x1 == i
#       #   for (j in n12) nijk[i, j] <- sum(d.x.i & d.y1 == 
#       #                                      j)
#       # }
#       t.log <- n * (nijk/tcrossprod(rowSums(nijk), colSums(nijk)))
#     }, {
#       dmS.1 <- dm[, S] + 1L
#       nijk <- array(0L, c(2, 2, 2))
#       for (i in n12) {
#         d.x.i <- d.x1 == i
#         for (j in n12) {
#           d.x.i.y.j <- d.x.i & d.y1 == j
#           for (k in n12) nijk[i, j, k] <- sum(d.x.i.y.j & 
#                                                 dmS.1 == k)
#         }
#       }
#       alt <- c(x, y, S)
#       c <- which(alt == S)
#       nik <- apply(nijk, c, rowSums)
#       njk <- apply(nijk, c, colSums)
#       nk <- colSums(njk)
#       t.log <- array(0, c(2, 2, 2))
#       if (c == 3) {
#         for (k in n12) t.log[, , k] <- nijk[, , k] * 
#             (nk[k]/tcrossprod(nik[, k], njk[, k]))
#       } else if (c == 1) {
#         for (k in n12) t.log[k, , ] <- nijk[k, , ] * 
#             (nk[k]/tcrossprod(nik[, k], njk[, k]))
#       } else {
#         for (k in n12) t.log[, k, ] <- nijk[, k, ] * 
#             (nk[k]/tcrossprod(nik[, k], njk[, k]))
#       }
#     }, {
#       dmS1.1 <- dm[, S[1]] + 1L
#       dmS2.1 <- dm[, S[2]] + 1L
#       nijk2 <- array(0L, c(2, 2, 2, 2))
#       for (i in n12) {
#         d.x.i <- d.x1 == i
#         for (j in n12) {
#           d.x.i.y.j <- d.x.i & d.y1 == j
#           for (k in n12) {
#             d.x.y.S1 <- d.x.i.y.j & dmS1.1 == k
#             for (l in n12) nijk2[i, j, k, l] <- sum(d.x.y.S1 & 
#                                                       dmS2.1 == l)
#           }
#         }
#       }
#       nijk <- array(0L, c(2, 2, 4))
#       for (i in n12) {
#         for (j in n12) nijk[, , 2 * (i - 1) + j] <- nijk2[, 
#                                                           , i, j]
#       }
#       nik <- apply(nijk, 3, rowSums)
#       njk <- apply(nijk, 3, colSums)
#       nk <- colSums(njk)
#       t.log <- array(0, c(2, 2, 4))
#       for (k in 1:4) t.log[, , k] <- nijk[, , k] * (nk[k]/tcrossprod(nik[, 
#                                                                          k], njk[, k]))
#     }, {
#       dmS1.1 <- dm[, S[1]] + 1L
#       dmS2.1 <- dm[, S[2]] + 1L
#       dmS3.1 <- dm[, S[3]] + 1L
#       nijk <- array(0L, c(2, 2, 8))
#       for (i1 in n12) {
#         d.x.i <- d.x1 == i1
#         for (i2 in n12) {
#           d.x.y.i12 <- d.x.i & d.y1 == i2
#           for (i3 in n12) {
#             d.x.y.S1 <- d.x.y.i12 & dmS1.1 == i3
#             for (i4 in n12) {
#               d.xy.S1S2 <- d.x.y.S1 & dmS2.1 == i4
#               for (i5 in n12) nijk[i1, i2, 4 * (i3 - 
#                                                   1) + 2 * (i4 - 1) + i5] <- sum(d.xy.S1S2 & 
#                                                                                    dmS3.1 == i5)
#             }
#           }
#         }
#       }
#       nik <- apply(nijk, 3, rowSums)
#       njk <- apply(nijk, 3, colSums)
#       nk <- colSums(njk)
#       t.log <- array(0, c(2, 2, 8))
#       for (k in 1:8) t.log[, , k] <- nijk[, , k] * (nk[k]/tcrossprod(nik[, 
#                                                                          k], njk[, k]))
#     }, {
#       dmS1.1 <- dm[, S[1]] + 1L
#       dmS2.1 <- dm[, S[2]] + 1L
#       dmS3.1 <- dm[, S[3]] + 1L
#       dmS4.1 <- dm[, S[4]] + 1L
#       nijk <- array(0L, c(2, 2, 16))
#       for (i1 in n12) {
#         d.x.i <- d.x1 == i1
#         for (i2 in n12) {
#           d.x.y.i12 <- d.x.i & d.y1 == i2
#           for (i3 in n12) {
#             d.x.y.S1 <- d.x.y.i12 & dmS1.1 == i3
#             for (i4 in n12) {
#               d.xy.S1S2 <- d.x.y.S1 & dmS2.1 == i4
#               for (i5 in n12) {
#                 d.xy.S1S2S3 <- d.xy.S1S2 & dmS3.1 == 
#                   i5
#                 for (i6 in n12) nijk[i1, i2, 8 * (i3 - 
#                                                     1) + 4 * (i4 - 1) + 2 * (i5 - 1) + 
#                                        i6] <- sum(d.xy.S1S2S3 & dmS4.1 == 
#                                                     i6)
#               }
#             }
#           }
#         }
#       }
#       nik <- apply(nijk, 3, rowSums)
#       njk <- apply(nijk, 3, colSums)
#       nk <- colSums(njk)
#       t.log <- array(0, c(2, 2, 16))
#       for (k in 1:16) t.log[, , k] <- nijk[, , k] * (nk[k]/tcrossprod(nik[, 
#                                                                           k], njk[, k]))
#     }, {
#       dmS1.1 <- dm[, S[1]] + 1L
#       dmS2.1 <- dm[, S[2]] + 1L
#       dmS3.1 <- dm[, S[3]] + 1L
#       dmS4.1 <- dm[, S[4]] + 1L
#       dmS5.1 <- dm[, S[5]] + 1L
#       nijk <- array(0L, c(2, 2, 32))
#       for (i1 in n12) {
#         d.x.i <- d.x1 == i1
#         for (i2 in n12) {
#           d.x.y.i12 <- d.x.i & d.y1 == i2
#           for (i3 in n12) {
#             d.x.y.S1 <- d.x.y.i12 & dmS1.1 == i3
#             for (i4 in n12) {
#               d.xy.S1S2 <- d.x.y.S1 & dmS2.1 == i4
#               for (i5 in n12) {
#                 d.xy.S1S2S3 <- d.xy.S1S2 & dmS3.1 == 
#                   i5
#                 for (i6 in n12) {
#                   d.xy.S1S2S3S4 <- d.xy.S1S2S3 & dmS4.1 == 
#                     i6
#                   for (i7 in n12) nijk[i1, i2, 16 * 
#                                          (i3 - 1) + 8 * (i4 - 1) + 4 * (i5 - 
#                                                                           1) + 2 * (i6 - 1) + i7] <- sum(d.xy.S1S2S3S4 & 
#                                                                                                            dmS5.1 == i7)
#                 }
#               }
#             }
#           }
#         }
#       }
#       nik <- apply(nijk, 3, rowSums)
#       njk <- apply(nijk, 3, colSums)
#       nk <- colSums(njk)
#       t.log <- array(0, c(2, 2, 32))
#       for (k in 1:32) t.log[, , k] <- nijk[, , k] * (nk[k]/tcrossprod(nik[, 
#                                                                           k], njk[, k]))
#     })
#   }
#   else {
#     nijk <- array(0L, c(2, 2, 1))
#     i <- d.x1[1]
#     j <- d.y1[1]
#     k <- NULL
#     lapply(as.list(S), function(x) {
#       k <<- cbind(k, dm[, x] + 1)
#       NULL
#     })
#     parents.count <- 1L
#     parents.val <- t(k[1, ])
#     nijk[i, j, parents.count] <- 1L
#     for (it.sample in 2:n) {
#       new.p <- TRUE
#       i <- d.x1[it.sample]
#       j <- d.y1[it.sample]
#       t.comp <- t(parents.val[1:parents.count, ]) == k[it.sample, 
#                                                        ]
#       dim(t.comp) <- c(lenS, parents.count)
#       for (it.parents in 1:parents.count) {
#         if (all(t.comp[, it.parents])) {
#           nijk[i, j, it.parents] <- nijk[i, j, it.parents] + 
#             1L
#           new.p <- FALSE
#           break
#         }
#       }
#       if (new.p) {
#         parents.count <- parents.count + 1L
#         if (verbose >= 2) 
#           cat(sprintf(" adding new parents (count = %d) at sample %d\n", 
#                       parents.count, it.sample))
#         parents.val <- rbind(parents.val, k[it.sample, 
#                                             ])
#         nijk <- abind(nijk, array(0, c(2, 2, 1)))
#         nijk[i, j, parents.count] <- 1L
#       }
#     }
#     nik <- apply(nijk, 3, rowSums)
#     njk <- apply(nijk, 3, colSums)
#     nk <- colSums(njk)
#     t.log <- array(0, c(2, 2, parents.count))
#     for (k in 1:parents.count) t.log[, , k] <- nijk[, , 
#                                                     k] * (nk[k]/tcrossprod(nik[, k], njk[, k]))
#   }
#   G2 <- sum(2 * nijk * log(t.log), na.rm = TRUE)
#   if (adaptDF && lenS > 0) {
#     zero.counts <- sum(nijk == 0L) + 4 * (2^lenS - dim(nijk)[3])
#     ndf <- max(1, df - zero.counts)
#     if (verbose) 
#       cat("adaptDF: (df=", df, ", zero.counts=", zero.counts, 
#           ") ==> new df = ", ndf, "\n", sep = "")
#     df <- ndf
#   }
#   pchisq(G2, df, lower.tail = FALSE)
# }
# 
# cpd_1 <- function(d.x1, d.y1){
#   nijk <- array(0L, c(2, 2))
#   for (i in n12) {
#     d.x.i <- d.x1 == i
#     for (j in n12) nijk[i, j] <- sum(d.x.i & d.y1 == 
#                                        j)
#   }
#   return(nijk)
# }
