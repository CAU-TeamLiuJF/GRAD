setwd("/public/home/liujf/workspace/xueyh/TempWork/grad")


options(warn = -1)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
library(rrBLUP)
library(lme4)
library(lme4qtl)

# pheno  Gmatrix
load("data/pig_data.rdata")
id.code <- read_delim("data/id.code.txt")

## gene expression
ge <- readRDS("data/yy_gene_fpkm_clean.RDS") %>% 
  rownames_to_column("id") %>% 
  left_join(id.code, by = "id") %>% 
  select(-id) %>% 
  filter(num %in% rownames(TT)) %>% 
  column_to_rownames("num") %>% 
  as.data.frame()

## filter 
gmat <- GG[rownames(TT), colnames(TT)]
ge <- ge[rownames(TT), ]

## genotype
geno <- fread("data/geno.raw", data.table = F) %>% 
  select(-c(1, 3:6), id = 2) 

pheno <- phe %>% rename(id = 1)

## dominance
library(AGHmatrix)

geno = geno %>% column_to_rownames("id") %>% as.matrix

GD = AGHmatrix::Gmatrix(geno, missingValue = -9, method = "Su")

## function to extract additive and dominance genetic effect in omics features
## x: gene expressions n*m 
## g1: additive Gmatrix
## g2: dominance Dmatrix
extract_genetic_eff <- function(x, g1, g2, pb = TRUE) {
    # 初始化结果存储
    temp_add <- matrix(NA, nrow = nrow(x), ncol = ncol(x)) %>% as.data.frame()
    temp_dom <- matrix(NA, nrow = nrow(x), ncol = ncol(x)) %>% as.data.frame()
    temp_res <- matrix(NA, nrow = nrow(x), ncol = ncol(x)) %>% as.data.frame()
    
    # 新增：存储方差组分的 data.frame
    var_components <- data.frame(
        trait = colnames(x)[-1],  # 假设第一列是ID
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
	
    for (i in 2:ncol(x)) {  # 假设第1列是ID
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
            # 提取随机效应值
            temp_add[same_id, i] <- ranef(tmp)[[1]][, 1]
            temp_dom[same_id, i] <- ranef(tmp)[[2]][, 1]
            temp_res[same_id, i] <- residuals(tmp)
            
            # 提取方差组分
            vc <- as.data.frame(VarCorr(tmp))
            add_var <- vc$vcov[1]       # 加性方差
            dom_var <- vc$vcov[2]      	# 显性方差
            res_var <- vc$vcov[3] 		# 残差方差
            total_var <- sum(add_var, dom_var, res_var)
            
            # 保存到方差组分表
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
        
        #if (i %% 100 == 0) cat("已完成: ", i, "/", ncol(x), "\n")
		
		if (pb) cli::cli_progress_update()
    }
    
    # 整理输出
    temp_add[, 1] <- temp_dom[, 1] <- temp_res[, 1] <- rownames(x)
    
	if (pb) cli::cli_process_done()
	
    return(list(
        add = temp_add,
        dom = temp_dom,
        res = temp_res,
        var_components = var_components
    ))
}

ge_id = ge %>% rownames_to_column("id")
rownames(ge_id) = ge_id$id

## calculate genetic effect in omics features
ge_eff = extract_genetic_eff(ge_id, GG, GD)
saveRDS(ge_eff, "data/yy_gene_fpkm_clean_eff.RDS")

## additive omics kernel
dt <- as.data.table(ge_eff[["add"]][,-1])
dt2 = dt[, which(colMeans(dt == 0, na.rm = TRUE) < 0.95), with = FALSE]
TT3 <- dt2 %>% 
  as.matrix() %>% 
  scale() %>% 
  { tcrossprod(.) / ncol(.) }
rownames(TT3)=colnames(TT3)=ge_eff[["add"]][,1]

## dominance omics kernel 
dt3 <- as.data.table(ge_eff[["dom"]][,-1])
dt4 = dt3[, which(colMeans(dt3 == 0, na.rm = TRUE) < 0.95), with = FALSE]  
TD <- dt4 %>% 
  as.matrix() %>% 
  scale() %>% 
  { tcrossprod(.) / ncol(.) }
rownames(TD)=colnames(TD)=ge_eff[["dom"]][,1]  

save(AA, GG, GD, TT, TT2, TT3, TD, pheno, geno, file = paste0("data/data_add_dom.rdata"))

## variance components
ge_id = ge %>% rownames_to_column("id")
rownames(ge_id) = ge_id$id
ge_eff = extract_genetic_eff(ge_id, GG, GD)
saveRDS(ge_eff, "data/yy_gene_fpkm_clean_eff_varcom.RDS")

