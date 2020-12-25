## Script to reproduce figure 5 of the manuscript:
## SPINT2 controls SARS-CoV-2 viral infection and is associated to disease severity

## Dependencies
library(magrittr)
library(tidyverse)
library(reshape2)
library(Seurat)
library(viridis)

markers <- c('TMPRSS2', 'ACE2', 
             'SPINT2')

############################################################################
## Expression of markers in Lung Adenocarcinoma

luad <- readRDS('~/sc/sars-cov2/data/laughney_processed_seu.rds')

luad$normal_vs_tumor <- factor(luad$normal_vs_tumor,
                               levels = c('TUMOR', 'NORMAL'))

pdf('/Users/carlosramirez/sc/sars-cov2/figures/luad_permissivity_genes.pdf')
DotPlot(subset(luad, abbv == 'EPITHELIAL'), 
        features = markers, 
        group.by = 'normal_vs_tumor',
        dot.scale = 40) +
        scale_color_viridis() +
        theme_bw(base_size = 9) +
        theme(panel.grid = element_blank(),
              axis.text = element_text(size = 12)) +
        xlab('') + ylab('') +
        coord_flip()
dev.off()


#########################################################################
## Colorectal Cancer sc-RNA seq evaluation
## Reference: https://www.nature.com/articles/ng.3818
epi <- readRDS('~/sc/sars-cov2/data/GSE81861_CRC/epi_seu.rds')

## Expression of comorbidity markers in normal vs tumor cells
## in Colorectal Cancer
pdf('~/sc/sars-cov2/figures/CRC_permissivity_genes_expression.pdf',
    width = 6, height = 8)
DotPlot(epi, features = markers, 
        group.by = 'normal_vs_tumor',
        dot.scale = 40) + 
        scale_color_viridis() +
        theme_bw(base_size = 9) +
        theme(text = element_text(color='black'),
              panel.grid = element_blank(),
              axis.text = element_text(size=18)) +
        coord_flip() +
        xlab('') + ylab('')
dev.off()


########################################################################
## Kidney Tumor sc-RNA Seq Behjati dataset
behjati <- readRDS('~/sc/sars-cov2/analysis/behjati2018_endothelial_seu.rds')


behjati$Main_category <- factor(behjati$Main_category,
                                levels = c('Tumor', 'Normal'))

pdf('~/sc/sars-cov2/figures/kidney_permissivity_genes_expression.pdf',
    width = 6, 
    height = 8)
DotPlot(behjati, 
        features = markers, 
        group.by = 'Main_category',
        dot.scale = 45) + 
        scale_color_viridis() +
        theme_bw(base_size = 9) +
        theme(text = element_text(color='black'),
              panel.grid = element_blank(),
              axis.text = element_text(size=18)) +
        coord_flip() +
        xlab('') + ylab('')
dev.off()

######################################################################
## Hepatocellular cancer

## Lu data on Hepatocellular carcinoma scRNA-Seq
url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE149614&format=file&file=GSE149614%5FHCC%2EscRNAseq%2ES71915%2Ecount%2Etxt%2Egz'

download.file(url = url, destfile = 'data/GSE149614_HCC.scRNAseq.S71915.count.txt.gz')


liver.mtx <- read.table('data/GSE149614_HCC.scRNAseq.S71915.count.txt.gz', 
                        sep = '\t', quote = '', header = TRUE)


## Definition of the seurat object
liver_seu <- CreateSeuratObject(
  counts = liver.mtx, 
  project = 'SARS-COV2', 
  assay = 'RNA', 
  min.cells = 1, 
  min.features = 1
)

## Seurat standard workflow
liver_seu <- SCTransform(liver_seu)
liver_seu <- FindVariableFeatures(liver_seu, nfeatures = 3000)
liver_seu <- RunPCA(liver_seu, features = VariableFeatures(liver_seu))
liver_seu <- RunUMAP(liver_seu, dims = 1:30)
liver_seu <- FindNeighbors(liver_seu, dims = 1:20)
liver_seu <- FindClusters(liver_seu, resolution = 3)

## Annotation of samples
liver_seu$'sample' <- gsub('_.*', '', colnames(liver_seu)) 
liver_seu$'selected_tumor_or_normal' <- grepl('N|T', liver_seu$sample)
liver_seu.ss <- subset(liver_seu, selected_tumor_or_normal == TRUE)
liver_seu.ss$'tumor_vs_normal' <- sapply(liver_seu.ss$sample, 
                                         function(x) ifelse(grepl('T', x),
                                                            'tumor', 'normal'))
liver_seu.ss <- FindNeighbors(liver_seu.ss, dims = 1:20)
liver_seu.ss <- FindClusters(liver_seu.ss, resolution = 1)

clusters <- c('0' = 'NK',
              '1' = 'NK',
              '2' = 'hepatocytes',
              '3' = 'NK',
              '4' = 'NK',
              '5' = 'kupffer cells',
              '6' = 'hepatocytes',
              '7' = 'endothelia',
              '8' = 'hepatocytes',
              '9' = 'kupffer cells',
              '10' = 'epithelia',
              '11' = 'NA',
              '12' = 'NA',
              '13' = 'NK',
              '14' = 'NA',
              '15' = 'NK',
              '16' = 'NA',
              '17' = 'endothelia',
              '18' = 'kupffer cells',
              '19' = 'hepatocytes',
              '20' = 'NK',
              '21' = 'NA',
              '22' = 'NK',
              '23' = 'NA',
              '24' = 'NA',
              '25' = 'NA',
              '26' = 'NA',
              '27' = 'hepatocytes',
              '28' = 'NA',
              '29' = 'NA',
              '30' = 'NA',
              '31' = 'NA',
              '32' = 'NA',
              '33' = 'hepatocytes',
              '34' = 'NK',
              '35' = 'NA',
              '36' = 'NA',
              '37' = 'NA',
              '38' = 'NA')

## Assigning clusters
liver_seu.ss$'cell_type' <- plyr::mapvalues(liver_seu.ss$seurat_clusters,
                                            from = names(clusters),
                                            to = clusters)

pdf('~/sc/sars-cov2/figures/figure_04_liver.pdf',
    width = 5, height = 5)
hepatocytes <- subset(liver_seu.ss, cell_type == 'hepatocytes')
DotPlot(hepatocytes, 
        features = c('SPINT2', 'ACE2', 'TMPRSS2'),
        group.by = 'tumor_vs_normal',
        dot.scale = 20) +
        scale_color_viridis() +
        theme(axis.text.x = element_text(angle = 45, 
                                         hjust = 1)) +
        xlab('') + ylab('') +
        coord_flip()
dev.off()


