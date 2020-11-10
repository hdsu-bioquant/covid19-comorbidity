### This script reproduces the supplementary figure 1

## Dependencies
library(dplyr)
library(Seurat)
library(ggplot2)

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

############################################################
## Visualization of the results
scores <- readRDS('analysis/hcl_scoring.rds')
head(scores)

my_theme <- theme_bw() +
                theme(panel.grid = element_blank(),
                      legend.position = 'none') 
        
## Cell type plot
scores_by_cell_type <- split(scores$permissivity_signature_score1, 
                                f = scores$celltype) 
mean_by_celltype <- lapply(scores_by_cell_type, mean)
order.celltype <- mean_by_celltype %>%
        unlist() %>%
        sort() %>%
        names

## Violin plot. Permissivity scores of cell types
pdf('figures/vlnplot_permissivity_score_bycelltype_hcl.pdf',
    width = 13)
scores %>%
        mutate(celltype = factor(celltype,
                                 levels = order.celltype)) %>%
        ggplot(aes(x=celltype, 
                   y=permissivity_signature_score1,
                   fill=celltype)) +
                geom_violin() +
                coord_flip() +
                my_theme +
                xlab('') + ylab('Permissivity signature score')
dev.off()

## Box plot
pdf('figures/boxplot_permissivity_score_bycelltype_hcl.pdf',
    width = 13)
scores %>%
        mutate(celltype = factor(celltype,
                                 levels = order.celltype)) %>%
        ggplot(aes(x=celltype, 
                   y=permissivity_signature_score1,
                   fill=celltype)) +
        geom_boxplot() +
        coord_flip() +
        my_theme +
        xlab('') + ylab('Permissivity signature score')
dev.off()

## By tissue type

## Getting order
scores <- mutate(scores, 
                 sample = gsub('Adult', '', sample))
scores_by_tissue <- split(scores$permissivity_signature_score1, 
                             f = scores$sample) 
mean_by_tissue <- lapply(scores_by_tissue, mean)
order.tissue <- mean_by_tissue %>%
        unlist() %>%
        sort() %>%
        names

## Violin plot. Permissivity scores by tissue
pdf('figures/vlnplot_permissivity_score_bytissue_hcl.pdf',
    width = 13)
scores %>%
        mutate(sample = factor(sample,
                                 levels = order.tissue)) %>%
        mutate(sample = gsub('Adult', '', sample)) %>%
        ggplot(aes(x=sample, 
                   y=permissivity_signature_score1,
                   fill=sample)) +
        geom_violin() +
        coord_flip() +
        my_theme +
        xlab('') + ylab('Permissivity signature score')
dev.off()

### Box plot. Permissivity scores by tissue
pdf('figures/boxplot_permissivity_score_by_tissue_hcl.pdf',
    width = 15)
scores %>%
        mutate(sample = factor(sample,
                               levels = order.tissue)) %>%
        ggplot(aes(x=sample, 
                   y=permissivity_signature_score1,
                   fill=sample)) +
        geom_boxplot() +
        coord_flip() +
        my_theme +
        xlab('') + ylab('Permissivity signature score')
dev.off()
