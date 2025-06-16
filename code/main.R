setwd("/public/home/liujf/workspace/xueyh/TempWork/grad/")
library(tidyverse)
library(rrBLUP)
library(lme4)
library(lme4qtl)

source("combine_matrix.R")
source("cross_validation.R")

## receive parameter
tr = as.numeric(commandArgs(t=T)[1])
cat("\n\n tr is : ", tr, "\n\n")

## load data
load("data/data_add_dom.rdata")

## expand the parameters
df = expand_grid(trait = c("age", "bf", "czs"), a = seq(0, 1, 0.05))
df = df[tr, ]

## print 
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


