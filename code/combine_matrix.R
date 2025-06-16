## single-step matrix for additive and dominance with Genotype and Gene expression data
## M1		n*n matrix
## M2       m*m matrix, n >= m
## a        number, weighting parameter
sshm <- function(M1, M2, a = 1) {
  id1 <- as.character(rownames(M1))
  id2 <- as.character(rownames(M2))
  id1 <- id1[!id1 %in% id2]
							    
  M1_11 <- M1[id1, id1]
  M1_12 <- M1[id1, id2]
  M1_21 <- M1[id2, id1]
  M1_22 <- M1[id2, id2]

  iM1_22 <- solve(M1_22)

  M2 <- M2[id2, id2]

  nind <- nrow(M2)
  avg_sum_1 <- sum(M2) / (nind * nind)
  avg_sum_2 <- sum(M1_22) / (nind * nind)
  avg_diag_1 <- sum(diag(M2)) / nind
  avg_diag_2 <- sum(diag(M1_22)) / nind
  sw <- (avg_sum_2 - avg_diag_2) / (avg_sum_1 - avg_diag_1)
  sm <- avg_diag_1 - sw * avg_diag_2
  M2 <- sw * M2 + sm

  M2 <- a * M2 + (1 - a) * M1_22

  H11 <- M1_11 + M1_12 %*% iM1_22 %*% (M2 - M1_22) %*% iM1_22 %*% M1_21
  H12 <- M1_12 %*% iM1_22 %*% M2
  H21 <- M2 %*% iM1_22 %*% M1_21
  H <- cbind(rbind(H11, H21), rbind(H12, M2))
  H <- H / mean(diag(H))

  return(H)
}