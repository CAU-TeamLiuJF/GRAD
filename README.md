# GRAD

## Model summary

We present a novel methods for integrating intermediate omics information into genomic prediction, which decomposes raw omics features into genetically regulated additive and dominance components while jointly modeling additive-dominance architectures at both genomic and intermediate omics levels.

## Requirements

R packages:

-`tidyverse`
 
-`lme4`
 
-`lme4qtl`
 
-`rrBLUP`

## Example data

There are required matrices for pig populations.

## How to run script

### Step 1: calculate the additive and dominance genomic-base matrix

### Step 2: calculate the additive and dominance omics-base matrix

> `R code/data_prepare.R` (Incorporated code in `step 1`)

### Step 3: prediction

> `R code/main.R`
