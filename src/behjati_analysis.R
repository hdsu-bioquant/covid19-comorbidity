## Dependencies
library(Seurat)
library(Matrix)
library(tidyverse)

############################################################################
## Auxiliary functions

## tweaked dot plot
tw_DotPlot <- function(seurat, 
                       features, 
                       dot_scale=10,
                       group_by){
        DotPlot(seurat, 
                features = features, 
                group.by = group_by, 
                dot.scale = dot_scale,
                cols = c('yellow', 'red')) + 
                theme(axis.text.x = element_text(angle = 90),
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank()) 
}

############################################################################
## Possible coreceptor of SARS-COV2

## interactome
interactome <- read.table('~/Downloads/interactome.tsv',
                          header = TRUE, 
                          stringsAsFactors = FALSE)
interactome <- interactome$PreyGene

## transcriptome
transcriptome <- read.table('~/Downloads/blanco_melo_transcriptome.DEG.tsv',
                            header = TRUE, 
                            stringsAsFactors = FALSE)
transcriptome <- transcriptome$GeneName

## Biogrid
biogrid <- read.table('~/Downloads/BIOGRID-CORONAVIRUS-3.5.184.tab3.txt',
                      header = TRUE, 
                      stringsAsFactors = FALSE, 
                      sep = '\t')
biogrid <- biogrid$Official.Symbol.Interactor.B

## Ashwin markers
ashwin <- c('ANPEP', 'DPP4', 'IFITM3', 
            'B2M', 'RPL36', 'PABPC1', 'RHOA')

## Union of all markers
markers <- unique(c(interactome, transcriptome, biogrid, ashwin,'ACE2', 'TMPRSS2'))


###############################################################################

## Reading Behjati data

## Count matrix
kidney <- readMM('~/Downloads/aat1699_DataS1/tableOfCounts.mtx')
## gene annotations
gene_anns <- read.table('~/Downloads/aat1699_DataS1/tableOfCounts_rowLabels.tsv',
                        header = T)
## cell annotations
anns <- readRDS('~/Downloads/cell_annot.RDS')
barcodes <- read.table('~/Downloads/aat1699_DataS1/tableOfCounts_colLabels.tsv',
                       header = TRUE)
anns <- data.frame(lapply(anns, as.character), stringsAsFactors = F)
anns <- cbind(anns, select(barcodes, DropletID, Barcode))

## annotation of dge matrix
colnames(kidney) <- anns$DropletID
rownames(kidney) <- gene_anns$Symbol
kidney <- kidney[!duplicated(rownames(kidney)), ]

################################################################################
## Filtering data

# exclude all the cells that did not pass QC ("Indistinct" were all the cells that did not pass QC according to annotation of the authors)
non_indistinct <- ! anns$Compartment == "Indistinct"

# subset cds (exclude Foetal cells), because we only want to compare Normal to Tumor cells
normal_vs_tumor <- with(anns, Main_category == 'Normal' | Main_category == 'Tumor' )

# subset cds to exclude immune cells and ureter samples
kidney_or_tumor <- sapply(anns$Category, 
                    function(x) x %in% c("Normal_mature_kidney", 
                                         "Kidney_tumour"))

## keeping only CCC
RCC <- grepl('RCC1|RCC2|VHL|pRCC', anns$Source)

## Dropping NA values
non_na <- ! is.na(anns$Main_category)

## performing subsetting
subset <- non_indistinct & normal_vs_tumor & kidney_or_tumor & non_na & RCC
sub_indexes <- (1:ncol(kidney))[subset]
kidney <- kidney[ ,  sub_indexes]

###########################################################################
## definition of the seurat object
behjati <- CreateSeuratObject(
        counts = kidney,
        min.cells = 1,
        min.features = 1,
        project = 'sars-cov2', 
        assay = 'RNA'
)


## Adding metadata
md <- anns[anns$DropletID %in% colnames(behjati), ]
behjati@meta.data <- cbind(behjati@meta.data, md)

## Normalization
behjati <- SCTransform(behjati)

## Saving seurat file
saveRDS(behjati, '~/sc/sars-cov2/analysis/behjati2018_seu.rds')

## Dimensional reduction
behjati <- RunPCA(behjati, npcs = 15, verbose = F)
behjati <- RunUMAP(behjati, reduction = 'pca', dims = 1:15)

## Explorative UMAPs
DimPlot(behjati, reduction = 'umap', group.by = 'Main_category')
DimPlot(behjati, reduction = 'umap', group.by = 'Cell_type2')
DimPlot(behjati, reduction = 'umap', group.by = 'Cell_type1')
DimPlot(behjati, reduction = 'umap', group.by = 'Category')
DimPlot(behjati, reduction = 'umap', group.by = 'Compartment')

## Inspecting normal vs tumor by cell type
behjati$'cell_type' <- paste(behjati$Category, 
                             behjati$Cell_type1, sep = ' + ')
DimPlot(behjati, 
        reduction = 'umap', 
        group.by = 'cell_type', 
        label = TRUE, repel = T) + NoLegend()


##########################################################################
## Normal vs tumor tissue DEG
Idents(behjati) <- behjati$Main_category
behjati.markers <- FindAllMarkers(behjati, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
behjati.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

## Extracting DEG genes
deg_mtx <- FetchData(behjati, 
                     vars = behjati.markers$gene) %>%
        as.matrix() %>% t
metadata <- behjati@meta.data

## Plotting heatmap
normal_vs_tumor_cols <- c('Normal' = 'green', 'Tumor' = 'red')
cell_type_cols <- c("ENDOTHELIAL" = 'salmon', "EPITHELIAL" = 'steelblue',
                    "FIBROBLAST" = 'purple', "PERICYTE" = 'aquamarine3',
                    "PMP" = 'chartreuse3')
col_anns <- HeatmapAnnotation(
        Sample_type = behjati$Main_category,
        col = list(Sample_type = normal_vs_tumor_cols))
H <- Heatmap(deg_mtx, 
             top_annotation = col_anns, 
             show_column_names = FALSE,
             column_title = 'Normal vs Tumor - Kidney DEG', 
             show_row_names = FALSE
)
H

############################################################################
## Fig 1. Intersection of DEG and sarscov-2 markers
deg_int_markers <- intersect(behjati.markers$gene, markers)


## Extracting intersection
deg_mtx <- FetchData(behjati, 
                     vars = deg_int_markers) %>%
        as.matrix() %>% t
metadata <- behjati@meta.data

## Defining cell types abbreviatures for labeling
behjati$'abbv' <- paste(behjati$Main_category, behjati$Cell_type1)
abbs <- c('Normal + Nep Others',
          'Normal + Nep Epi',
          'Normal + Endo',
          'Tumor + Nor',
          'RCC',
          'Tumor Endo',
          'Wilms + Fibro',
          'pRCC', 
          'Wilms Tumor')
behjati$abbv <- plyr::mapvalues(behjati$abbv,
                                from = unique(behjati$abbv),
                                to = abbs)
## Plotting heatmap
normal_vs_tumor_cols <- c('Normal' = 'green', 'Tumor' = 'red')
cell_type_cols <- c('salmon', 'steelblue', 'purple', 
                    'aquamarine3', 'chartreuse3', 'chocolate3',
                    'cornsilk2', 'brown3', 'darkgoldenrod1')
names(cell_type_cols) <- abbs
col_anns <- HeatmapAnnotation(
        Sample_type = behjati$Main_category,
        cell_type = behjati$abbv, 
        col = list(Sample_type = normal_vs_tumor_cols),
                   cell_type = cell_type_cols)
H <- Heatmap(deg_mtx, 
             top_annotation = col_anns, 
             show_column_names = FALSE,
             column_title = 'Normal vs Tumor - Kidney intersect(DEG, SARSCOV2 markers)'
)
H

writeLines(deg_int_markers, 
           '~/sc/sars-cov2/analysis/deg_int_markers_normal_vs_tumor_kidney.txt')


################################################################################
## Fig 2. DEG expressed genes dot plot 

## Reproducing Max results
subset(behjati, abbv != 'pRCC') %>%
        subset(abbv != 'Wilms Tumor') %>% 
          tw_DotPlot(
             features = c(deg_int_markers, 'ACTB'), 
             group_by = 'abbv')

subset(behjati, abbv != 'pRCC') %>%
        subset(abbv != 'Wilms Tumor') %>% 
      tw_DotPlot( 
           features = c('ACE2', 'TMPRSS2', 'NEU1', 'NPTX1'), 
           group_by = 'cell_type')


subset(behjati, abbv != 'pRCC') %>%
        subset(abbv != 'Wilms Tumor') %>% 
        tw_DotPlot( 
                features = deg_int_markers, 
                group_by = 'abbv')

#################################################################################
## Figure 3. DEG by cell type
Idents(behjati) <- behjati$abbv
behjati.markers <- FindAllMarkers(behjati, 
                                  only.pos = TRUE, 
                                  min.pct = 0.25, 
                                  logfc.threshold = 0.25)
behjati.markers.top <- behjati.markers %>% 
                           group_by(cluster) %>% 
                           top_n(n = 20, wt = avg_logFC)

behjati.markers.sp <- split(x = behjati.markers.top, 
                            f = behjati.markers.top$cluster)
deg_int_markers_ct <- lapply(behjati.markers.sp, 
                            function(x) intersect(x$gene, markers)) %>%
                              unlist()

kp <- subset(behjati, abbv != 'pRCC') %>%
        subset(abbv != 'Wilms Tumor') %>% 
        tw_DotPlot( 
                features = c(deg_int_markers_ct, 
                             c('ACTB','ACE2', 'TMPRSS2')), 
                group_by = 'abbv')

kp$data$id <- factor(kp$data$id, 
                     levels = c('Normal + Nep Others',
                                'Normal + Nep Epi', 
                                'Normal + Endo',
                                'RCC',
                                'Tumor + Nor',
                                'Tumor Endo',
                                'Wilms + Fibro'))
kp
