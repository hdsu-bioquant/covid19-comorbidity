## This script reproduces the supplementary figure 1

## Dependencies
library(dplyr)
library(Seurat)

## Initial settings
project_dir <- '/Users/carlosramirez/sc/diabetes_COVID19_comorbidity/'
setwd(project_dir)

##################################################################
##                                                              ##
##      Scoring permissivity signature in HCL cells             ##
##                                                              ##
##################################################################

## Reading HCL data
hcl_path <- '~/sc/serine_proteases_regulation/data/hcl/hcl_normalized_adult_seu.rds'
hcl <- readRDS(hcl_path)

## Subsetting to adult tissue
hcl$'adult' <- grepl('adult', tolower(hcl$orig.ident)) &
        ! ( grepl('fetal', tolower(hcl$celltype)))
hcl <- subset(hcl, adult == TRUE)

## Reading permissivity signature
perm_signature <- read.table(
        '~/sc/sars-cov2/data/permissivity_signature_annotated.tsv',
        header = TRUE
)
head(perm_signature)

## Cell scoring 
hcl <- AddModuleScore(
        hcl, 
        features = list(perm_signature=perm_signature), 
        nbin = 100
)
