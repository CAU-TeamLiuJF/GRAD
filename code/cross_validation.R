## cross validation
## pheno 	n*m data.frame,id in first colunmn
## pid 		trait name or colunmn number for analysised trait
## G 		n*n matrix for first random effect
## G2 		n*n matrix, optional, for second random effect
## fold		number, folds for CV 
## rep      number, repeat times for CV
## pb		bool, whether show the progress bar
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
  require(lme4)
  require(lme4qtl)
  
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
		require(rrBLUP)
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