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

###################################################################
##                                                               ##
##      SPINT2 and TMPRSS2 expression in Calu-3 and H1299        ##
##      cells                                                    ##
##                                                               ##
###################################################################

## Reading data
calu <- readRDS('../sars-cov2/data/200408.Seurat_Calu_CoV_1000_Merged.rds')
calu$strain[calu$strain == 'nan'] <- 'mock'
calu.nin <- subset(calu, orig.ident %in% c("Calu3-mock-4h-A", "Calu3-mock-4h-B") )
rm(calu)
h1299 <- readRDS('../sars-cov2/data/200406.Seurat_H1299_CoV_1000_Merged.rds')
h1299.nin <- subset(h1299, orig.ident %in% c("H1299-mock-4h-A", "H1299-mock-4h-B"))
rm(h1299)
calu.h1299.nin <- merge(calu.nin, h1299.nin)
rm(calu.nin, h1299.nin)

## Normalization and dimension reduction
calu.h1299.nin <- SCTransform(calu.h1299.nin)
calu.h1299.nin <- FindVariableFeatures(calu.h1299.nin, nfeatures = 3000)
calu.h1299.nin <- RunPCA(calu.h1299.nin, 
                         features = VariableFeatures(calu.h1299.nin))
ElbowPlot(calu.h1299.nin)
calu.h1299.nin <- RunUMAP(calu.h1299.nin, reduction = 'pca', dims = 1:10)

calu.h1299.nin$'cell_line' <- sapply(
    calu.h1299.nin$orig.ident, function(x)
    ifelse(grepl('Calu', x), 'Calu-3', 'H1299')
)
calu.h1299.nin$cell_line %>% table()

DotPlot(calu.h1299.nin,
        group.by = 'cell_line', 
        features = c('TMPRSS2', 'ACE2', 'SPINT2', '')
)
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
        sort(decreasing = TRUE) %>%
        names

## selecting top n cell types
n <- 30

## Violin plot. Permissivity scores of cell types
pdf('figures/vlnplot_permissivity_score_bycelltype_hcl.pdf',
    width = 13)
scores %>%
        filter(celltype %in% order.celltype[1:n]) %>%
        mutate(celltype = factor(celltype,
                                 levels = order.celltype[1:n])) %>%
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
    width = 13, height = 10)
scores %>%
        filter(celltype %in% order.celltype[1:n]) %>%
        mutate(celltype = factor(celltype,
                                 levels = order.celltype[1:n][n:1])
               ) %>%
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

############################################################
##                                                        ##
##      Negative correlation of SPINT2 and viral          ##
##      quantifications in Calu-3 and Caco-2              ##
##                                                        ##
############################################################

## Calu-3

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
## Caco proteomic data

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
spint2_viral_cor_traslatome <- traslatome %>%
        ggplot(aes(x=SPINT2, y=P0DTC2)) +
                geom_point() +
                geom_smooth(method = 'lm') +
                my_theme +
                ylab('P0DTC2 | Spike glycoprotein SARS-CoV2')
spint2_viral_cor_traslatome
dev.off()

pdf('figures/spint2_viral_correlation_cell_lines.pdf',
    width = 12)
gridExtra::grid.arrange(
        spint2_viral_expression_cor,
        spint2_viral_cor_traslatome,
        ncol = 2
)
dev.off()

##########################################################
## CDX2 binds to SPINT2 and TMPRSS2 during the embryonic
## development

markers <- c('SPINT2', 'TMPRSS2', 'SPINT1',
             'ST14')

chipseq <- read_xls('data/41598_2017_16009_MOESM3_ESM.xls',
                    sheet = 'Cdx2 Association Score', 
                    skip = 1)

pdf('figures/cdx2_chipseq_mouse_embryo.pdf')
chipseq %>%
    arrange(desc(`Cdx2 Score`)) %>%
    mutate(rank=1:nrow(chipseq)) %>%
    mutate(highlight=ifelse(rank < 10 | 
                                toupper(Symbol) %in% markers,
                            TRUE, FALSE)) %>%
    mutate(gene_label = ifelse(highlight == TRUE,
                               Symbol, '')) %>%
    ggplot(aes(x = rank,
               y = `Cdx2 Score`,
               colour = highlight,
               label = gene_label)) +
            geom_point() +
            scale_color_manual(values=c('black',
                                        'red')) +
            geom_text_repel() +
            theme_bw() +
            theme(panel.grid = element_blank(),
                  legend.position = 'none') 
dev.off()

