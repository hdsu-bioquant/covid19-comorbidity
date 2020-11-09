### This script reproduces the supplementary figure 1
#â
## Dependencies
library(dplyr)
library(Seurat)

## Initial settings
project_dir <- '/media/ag-cherrmann/cramirez/covid19-comorbidity/'
setwd(project_dir)

##################################################################
##                                                              ##
##      Scoring permissivity signature in HCL cells             ##
##                                                              ##
##################################################################

## Reading HCL data
hcl_path <- 'data/hcl_normalized_adult_seu.rds'
hcl <- readRDS(hcl_path)

## Subsetting to adult tissue
hcl$'adult' <- grepl('adult', tolower(hcl$orig.ident)) &
        ! ( grepl('fetal', tolower(hcl$celltype)))
hcl <- subset(hcl, adult == TRUE)
hcl

## Reading permissivity signature
perm_signature <- read.table(
        'analysis/permissivity_signature_annotated.tsv',
        header = TRUE
)
head(perm_signature)
indxs <- 0.25*nrow(perm_signature) %>% as.integer()
indxs
signature <- perm_signature$gene[1:indxs]

## Cell scoring 
hcl <- AddModuleScore(
        hcl, 
        features = list(signature),
        name = 'permissivity_signature_score', 
        nbin = 100
)
#mean(hcl$'permissivity_signature_score')

## Saving results
saveRDS(hcl@meta.data, 'analysis/hcl_scoring.rds')
write.table(hcl@meta.data, 
            'analysis/hcl_scoring.tsv', 
            sep='\t')
