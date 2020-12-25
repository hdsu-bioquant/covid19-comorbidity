### This script reproduces the supplementary figure 1

## Dependencies
library(dplyr)
library(Seurat)
library(ggplot2)
library(readxl)
library(ggrepel)

## Initial settings
project_dir <- '/Users/carlosramirez/sc/covid19-comorbidity/'
#project_dir <- '/media/ag-cherrmann/cramirez/covid19-comorbidity/'
setwd(project_dir)

######################################################################
##                                                                  ##
##      Supp Figure 1A. Negative correlation of SPINT2 to viral     ##
##      expression in Calu-3 cell lines                             ##
##                                                                  ##
######################################################################

## Linear correlation to viral gene expression
get_sign_scores <- function(seurat_object,
                            gene_set){
        gene.subset <- FetchData(seurat_object, gene_set)
        signature.scores <- rowSums(gene.subset)
        return(signature.scores)       
}

## Reading Calu data
calu <- readRDS('~/sc/sars-cov2/data/200408.Seurat_Calu_CoV_1000_Merged.rds')

## Subsetting infected at latter times post-infection (24)
calu.scov2.inf <- subset(calu, strain == 'SARSCoV2' & 
                                 infect == 'infected' &
                                 orig.ident %in% c("Calu3-S2-12h-A", 
                                                   "Calu3-S2-12h-B")
)
rm(calu)

## sars-cov genes 
scov2_genes <- grep('CoV2', rownames(calu.scov2.inf), value = T)

## normalization
calu.scov2.inf <- NormalizeData(calu.scov2.inf) %>% ScaleData()

## sars-cov2 expression
calu.scov2.inf$'sarscov2' <- get_sign_scores(calu.scov2.inf,
                                             gene_set = scov2_genes)

ge.df <- FetchData(calu.scov2.inf,
                   vars = c('sarscov2', 'SPINT2', 'TMPRSS2', 'ACE2'),slot = )

jpeg('figures/spint2_viral_expression_cor.jpeg')
spint2_viral_expression_cor <- ge.df %>% 
        ggplot(aes(x=SPINT2, y=sarscov2)) +
        geom_point(alpha=0.5) +
        geom_smooth(method = 'lm') +
        theme_classic() +
        ylab('Cumulative Sum of Viral Log2 Normalized Genes')
spint2_viral_expression_cor
dev.off()

################################################################
## Caco-2 translation rates

dir.create('data/')
dir.create('data/bojkova2020')

## Download data
file_url <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2332-7/MediaObjects/41586_2020_2332_MOESM2_ESM.xlsx'
file_path <- 'data/bojkova2020/41586_2020_2332_MOESM2_ESM.xlsx'
if ( ! file.exists(file_path) ) {
        download.file(url = file_url,
                      destfile = file_path)
}

## SARS-COV2 infected caco cell lines data
caco.traslatome <- read_xlsx(file_path, sheet = 1)

dim(caco.traslatome)

## Correlation to SARS-COV2 genes en caco traslatome
names(caco.traslatome)
## Adding gene names to viral proteins
caco.traslatome$'gene_name' <- sapply(1:nrow(caco.traslatome),
                                      function(i) ifelse(caco.traslatome$`Species Names01`[i] == 'Wuhan seafood market pneumonia virus OX=2697049',
                                                         caco.traslatome$Accession[i],
                                                         caco.traslatome$`Gene Symbol01`[i])
)
caco.traslatome <- caco.traslatome[!duplicated(caco.traslatome$gene_name),]
caco.traslatome <- caco.traslatome[!is.na(caco.traslatome$gene_name), ]
caco.traslatome <- as.data.frame(caco.traslatome)
rownames(caco.traslatome) <- caco.traslatome$gene_name

## Processing table
traslatome <- caco.traslatome[, grep('Virus', names(caco.traslatome))] %>%
                        t() %>% as.data.frame()

## Visualization
pdf('figures/spint2_viral_correlation_traslatome.pdf')
traslatome %>%
        ggplot(aes(x=SPINT2, y=P0DTC2)) +
                geom_point() +
                geom_smooth(method = 'lm') +
                my_theme +
                ylab('P0DTC2 | Spike glycoprotein SARS-CoV2')
dev.off()
