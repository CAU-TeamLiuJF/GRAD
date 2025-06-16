## function to extract additive and dominance genetic effect in omics features
## x: gene expressions n*m 
## g1: additive Gmatrix
## g2: dominance Dmatrix
## pb: show progress bar
extract_genetic_eff <- function(x, g1, g2, pb = TRUE) {
	require(cli)
	require(lme4)
	require(lme4qtl)
	
    # init data structure
    temp_add <- matrix(NA, nrow = nrow(x), ncol = ncol(x)) %>% as.data.frame()
    temp_dom <- matrix(NA, nrow = nrow(x), ncol = ncol(x)) %>% as.data.frame()
    temp_res <- matrix(NA, nrow = nrow(x), ncol = ncol(x)) %>% as.data.frame()
    
    var_components <- data.frame(
        trait = colnames(x)[-1], 
        add_var = NA,
        dom_var = NA,
        res_var = NA,
        add_prop = NA,
        dom_prop = NA,
        res_prop = NA
    )
    
    rownames(temp_add) <- rownames(temp_dom) <- rownames(temp_res) <- rownames(x)
    
    same_id = intersect(x[, "id"], colnames(g1))
    same_id = intersect(same_id, colnames(g2))
    x = x[same_id, ]
    G = g1[same_id, same_id]
    D = g2[same_id, same_id]
    
	if (pb) cli::cli_progress_bar("Progress", total = (ncol(x) - 1) )
	
    for (i in 2:ncol(x)) { 
        pheno = x[, c(1, 1, i)] %>% 
            rename(id = 1, id2 = 2, y = 3)
        
        tmp <- tryCatch(
            relmatLmer(
                y ~ (1 | id) + (1 | id2), 
                data = pheno, 
                relmat = list(id = G, id2 = D)
            ),
            error = function(e) {
                message("\nColunme ", i, " error: ", e$message)
                return(NULL)
            }
        )
        
        if (!is.null(tmp)) {
          
            temp_add[same_id, i] <- ranef(tmp)[[1]][, 1]
            temp_dom[same_id, i] <- ranef(tmp)[[2]][, 1]
            temp_res[same_id, i] <- residuals(tmp)
                        
            vc <- as.data.frame(VarCorr(tmp))
            add_var <- vc$vcov[1]       
            dom_var <- vc$vcov[2]      	
            res_var <- vc$vcov[3] 		
            total_var <- sum(add_var, dom_var, res_var)
            
            var_components[i-1, ] <- c(
                colnames(x)[i], 
                add_var, dom_var, res_var,
                add_var / total_var,
                dom_var / total_var,
                res_var / total_var
            )
        } else {
            var_components[i-1, ] <- c(colnames(x)[i], rep(NA, 6))
        }
		
		if (pb) cli::cli_progress_update()
    }
    
    # output
    temp_add[, 1] <- temp_dom[, 1] <- temp_res[, 1] <- rownames(x)
    
	if (pb) cli::cli_process_done()
	
    return(list(
        add = temp_add,
        dom = temp_dom,
        res = temp_res,
        var_components = var_components
    ))
}