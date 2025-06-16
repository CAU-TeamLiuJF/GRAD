setwd("/public/home/liujf/workspace/xueyh/TempWork/grad")
options(warn = -1)
library(tidyverse)
library(data.table)
library(rrBLUP)
library(lme4)
library(lme4qtl)
library(AGHmatrix)

source("extract_genetic_effect.R")


## read data
pheno <- read_delim("data/id.code.txt") %>% rename(id = 1)

## gene expression
ge <- readRDS("data/yy_gene_fpkm_clean.RDS") %>% 
  rownames_to_column("id") %>% 
  left_join(id.code, by = "id") %>% 
  select(-id) %>% 
  filter(num %in% rownames(TT)) %>% 
  column_to_rownames("num") %>% 
  as.data.frame()

## genotype
geno <- fread("data/geno.raw", data.table = F) %>% 
  select(-c(1, 3:6), id = 2) 


## filter for same ids 
gmat <- GG[rownames(TT), colnames(TT)]
ge <- ge[rownames(TT), ]


## kinship matrix
geno = geno %>% column_to_rownames("id") %>% as.matrix

GA = AGHmatrix::Gmatrix(geno, missingValue = -9, method = "VanRaden")
GD = AGHmatrix::Gmatrix(geno, missingValue = -9, method = "Su")

## extract genetic effect in omics features
ge_id = ge %>% rownames_to_column("id")
rownames(ge_id) = ge_id$id
ge_eff = extract_genetic_eff(ge_id, GG, GD)

## additive omics kernel
dt <- as.data.table(ge_eff[["add"]][,-1])
dt2 = dt[, which(colMeans(dt == 0, na.rm = TRUE) < 0.95), with = FALSE]
TA <- dt2 %>% 
  as.matrix() %>% 
  scale() %>% 
  { tcrossprod(.) / ncol(.) }
rownames(TA)=colnames(TA)=ge_eff[["add"]][,1]

## dominance omics kernel 
dt3 <- as.data.table(ge_eff[["dom"]][,-1])
dt4 = dt3[, which(colMeans(dt3 == 0, na.rm = TRUE) < 0.95), with = FALSE]  
TD <- dt4 %>% 
  as.matrix() %>% 
  scale() %>% 
  { tcrossprod(.) / ncol(.) }
rownames(TD)=colnames(TD)=ge_eff[["dom"]][,1]  

save(AA, GG, GD, TT, TA, TD, pheno, file = paste0("data/data_add_dom.rdata"))

