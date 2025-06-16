.libPaths("/public/home/liujf/workspace/xueyh/software/r4.3/lib/R/library")

setwd("/public/home/liujf/workspace/xueyh/TempWork/20250326_pig_transcriptome/")
library(tidyverse)
library(rrBLUP)
library(lme4)
library(lme4qtl)

# receive parameter
tr = as.numeric(commandArgs(t=T)[1])
cat("\n\n tr is : ", tr, "\n\n")

load("data/data_add_dom.rdata")

## single-step matrix for additive and dominance with Genotype and Gene expression data
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


## cross validation
crva <- function(
  pheno,
  pcol = 2,
  pid = 1,
  G,
  G2 = NULL,
  fold = 5,
  rep = 2,
  pb = FALSE,
  ...
) {
  if (is.null(pheno)) return(NULL)
  #G=G_sig;G2=G_res
  phenos <- pheno %>%
    dplyr::select(pid = pid, pcol) %>%
    filter(pid %in% rownames(G)) %>%
    mutate(gid = pid, gid = as.character(gid), gid2 = gid) %>%
    as.data.frame()

  G <- G[phenos$gid, phenos$gid]
  if (!is.null(G2)) G2 <- G2[phenos$gid, phenos$gid]

  results <- data.frame()
  namess <- ifelse(is.numeric(pcol), colnames(phenos)[pcol], pcol)
  if (pb) cli::cli_progress_bar(namess, total = fold * rep)
  for (i in seq_len(rep)) {
    for (j in seq_len(fold)) {
      set.seed(i)
      phenos$partition <- sample(
        seq_len(fold),
        size = nrow(phenos),
        replace = TRUE,
        prob = c(rep((1 / fold), times = fold))
      )
      phenos$yNA <- phenos[, pcol]
      nas <- phenos$partition == j
      phenos$yNA[nas] <- NA

      if (!is.null(G2)) {
        goo <- G[!nas, !nas]
        gno <- G[nas, !nas]
        too <- G2[!nas, !nas]
        tno <- G2[nas, !nas]

        res_gt <- relmatLmer(
          yNA ~ (1 | gid) + (1 | gid2),
          data = phenos[!nas, ],
          relmat = list(gid = G, gid2 = G2)
        )
        u1 <- as.matrix(
          t(res_gt@optinfo$relmat$relfac$gid) %*% as.matrix(ranef(res_gt)$gid)
        )[rownames(goo), 1]
        u2 <- as.matrix(
          t(res_gt@optinfo$relmat$relfac$gid2) %*% as.matrix(ranef(res_gt)$gid2)
        )[rownames(too), 1]

        pred_lme4qtl_gt <- gno %*%
          MASS::ginv(goo) %*%
          u1 +
          tno %*% MASS::ginv(too) %*% u2
      } else if (is.null(G2)) {
        ## rrBLUP - G
        res_k <- mixed.solve(phenos$yNA, K = G)

        pred_lme4qtl_gt <- res_k$u[nas]
      }

      # save result
      results <- rbind(
        results,
        data.frame(
          rep = i,
          fold = j,
          cor = cor(
            phenos[nas, pcol],
            pred_lme4qtl_gt,
            use = 'pairwise.complete.obs'
          ),
          bias = lm(phenos[nas, pcol] ~ pred_lme4qtl_gt)$coefficients[2]
        )
      )
      if (pb) cli::cli_progress_update()
    }
  }

  if (pb) cli::cli_process_done()
  return(results)
}

df = expand_grid(trait = c("age", "bf", "czs"), a = seq(0, 1, 0.05))
df = df[tr, ]

cat("\n tr is ", tr, "\n")
cat(" trait is ", df$trait[1], "\n")
cat(" a is ", df$a[1], "\n")


## mssBLUP
HA = sshm(GG, TT, df$a[1])

ssh = crva(
  pheno = pheno,
  pcol = df$trait[1],
  G = HA,
  fold = 5,
  rep = 10,
  pb = T
) %>%
  mutate(trait = df$trait[1])
  
if (!is.null(ssh)) {write_csv(ssh, paste0(getwd(), "/result/ss/ssh_", tr, ".csv"))}


## GRA
HA = sshm(GG, TT3, df$a[1])

ssha2 = crva(
  pheno = pheno,
  pcol = df$trait[1],
  G = HA,
  fold = 5,
  rep = 10,
  pb = T
) %>%
  mutate(trait = df$trait[1])

if (!is.null(ssha2)) {write_csv(ssha2, paste0(getwd(), "/result/ss_add2/ssha2_", tr, ".csv"))}


## GRAD
### with optimal weights in model GRA
if(df$trait == "bf") {HA = sshm(GG, TT2, 0)} else {HA = sshm(GG, TT2, 1)}
HD = sshm(GD, TD, df$a[1])

sshd = crva(
  pheno = pheno,
  pcol = df$trait[1],
  G = HA,
  G2 = HD,
  fold = 5,
  rep = 10,
  pb = T
) %>%
  mutate(trait = df$trait[1])


if (!is.null(sshd)) {write_csv(sshd, paste0(getwd(), "/result/ss_add_dom2/grad_", tr, ".csv"))}


## GBLUP
g_age = crva(pheno = pheno, pcol = "age", G = GG, fold = 5, rep = 10, pb = T)
g_bf = crva(pheno = pheno, pcol = "bf", G = GG, fold = 5, rep = 10, pb = T)
g_czs = crva(pheno = pheno, pcol = "czs", G = GG, fold = 5, rep = 10, pb = T)
write_csv(g_age, "result//gblup//age_gblup.csv")
write_csv(g_bf, "result//gblup//bf_gblup.csv")
write_csv(g_czs, "result//gblup//czs_gblup.csv")


## GDBLUP
gd_age = crva(pheno = pheno, pcol = "age", G = GG, G2 = GD, fold = 5, rep = 10, pb = T)
gd_bf = crva(pheno = pheno, pcol = "bf", G = GG, G2 = GD, fold = 5, rep = 10, pb = T)
gd_czs = crva(pheno = pheno, pcol = "czs", G = GG, G2 = GD, fold = 5, rep = 10, pb = T)
write_csv(gd_age, "result//gblup//age_gdblup.csv")
write_csv(gd_bf, "result//gblup//bf_gdblup.csv")
write_csv(gd_czs, "result//gblup//czs_gdblup.csv")


