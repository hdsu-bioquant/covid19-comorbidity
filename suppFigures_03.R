## Script to reproduce supplementary figure 3 of the manuscript
## SPINT2 controls SARS-CoV-2 viral infection and is associated to disease severity

## Dependencies
library(magrittr)
library(tidyverse)
library(reshape2)
library(Seurat)
library(viridis)

markers <- c('TMPRSS2', 'ACE2', 
             'SPINT2')

#######################################################################
##                                                                   ## 
##      Supp figure 3A. SPINT2 gene expression in tumor vs normal    ##
##      samples from TCGA                                            ##
##                                                                   ##
#######################################################################

## Reading preprocessed TCGA TPM data
tcga <- readRDS( 
        '/Users/carlosramirez/sc/sars-cov2/data/tcga_rpm/tcga_norm_vs_tumor_matched.rds'
)

## Reading gene annotations
gene_anns <- read.table(
        '/Users/carlosramirez/sc/sars-cov2/data/tcga_rpm/gencode.v23.annotation.gene.probemap',
        sep = '\t', 
        stringsAsFactors = FALSE, 
        header = TRUE
)
head(gene_anns)

## SPINT2 expression
marker <- 'SPINT2'
plots <- lapply(tcga, function(x) {
        x %>% as.data.frame() %>%
                mutate(id=rownames(x)) %>%
                merge(gene_anns) %>%
                filter(gene == marker) %>%
                select(matches('Normal|Cancer')) %>%
                t() %>%
                set_colnames('value') %>%
                as.data.frame() %>%
                rownames_to_column(var='tissue_type') %>%
                mutate(sample=gsub('-Cancer|-Normal', '', tissue_type),
                       tissue_type = ifelse(grepl('Cancer', tissue_type), 
                                            'cancer', 'normal')) %>%
                ggboxplot(x = "tissue_type", y = "value",
                          color = "tissue_type",
                          palette = "jco") +
                stat_compare_means(label.x = 1.25, 
                                   comparisons = list(c('cancer', 'normal'))) +
                theme(legend.position = 'none') +
                xlab('') + ylab(paste0(marker, ' gene expression TPM')) +
                ylim(0,14)
})

plots <- lapply(1:length(plots), function(i) plots[[i]] + ggtitle(label = names(plots)[i]) )
names(plots) <- names(tcga)
plots <- plots[names(plots) != 'OV']

pdf('/Users/carlosramirez/sc/sars-cov2/figures/TPM_TCGA_spint2.pdf',
    width = 10, height = 5)
gridExtra::grid.arrange(grobs = plots, ncol=5)
dev.off()

#######################################################################
##                                                                   ##
##      Supp figure 3B. Profiling of Hepatocellular carcinoma        ##
##      cells from Lu Y et al, 2020                                  ##
##                                                                   ##
#######################################################################

url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE149614&format=file&file=GSE149614%5FHCC%2EscRNAseq%2ES71915%2Ecount%2Etxt%2Egz'

download.file(url = url, destfile = 'data/GSE149614_HCC.scRNAseq.S71915.count.txt.gz')

liver.mtx <- read.table('data/GSE149614_HCC.scRNAseq.S71915.count.txt.gz', 
                        sep = '\t', quote = '', header = TRUE)

liver.mtx  <- readRDS('data/GSE149614_HCC.scRNAseq.S71915.count.rds')

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
DimPlot(liver_seu.ss, group.by = 'sample')
DimPlot(liver_seu.ss, group.by = 'tumor_vs_normal') 
DimPlot(liver_seu.ss, group.by = 'seurat_clusters')     

DotPlot(liver_seu.ss, 
        features = markers, 
        group.by = 'seurat_clusters') +
        scale_color_viridis() +
        theme(axis.text.x = element_text(angle = 45, 
                                         hjust = 1)) +
        xlab('') + ylab('')

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

DimPlot(liver_seu.ss, 
        group.by = 'cell_type') 

#######################################################################
##                                                                   ##
##      Supp figure 3C. SPINT2 gene expression in pancreatic         ##
##      cells from diabetic or control patients using                ##
##      data from Segerstolpe A. et al, 2016 scRNA-Seq               ##
##                                                                   ##
#######################################################################

## Data preprocessing
pancreas <- readRDS("~/sc/sars-cov2/data/E-MTAB-5061_seuobj_rpkm_qc.rds")
pancreas <- SCTransform(pancreas, 
                        vars.to.regress = "nCount_RNA", 
                        verbose = FALSE)

pancreas_alpha <- subset(pancreas, subset = cell_type == "alpha")
pancreas_alpha$condition <- factor(pancreas_alpha$condition, 
                                   levels = c('healthy_control', 'disease'))

pdf('figures/supp_figure_04_a_diabetes.pdf')
VlnPlot(pancreas, 
        features = 'SPINT2', 
        split.by = 'condition', 
        group.by = 'cell_type')
dev.off()





