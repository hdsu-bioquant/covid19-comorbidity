## This script contains the pre-processing of the data
## in the manuscript

## Dependencies
library(Seurat)
library(ggplot2)
library(tidyverse)
library(rhdf5)
library(Matrix)

setwd('/Users/carlosramirez/sc/sars-cov2/')
dir.create('data/')

###################################################################
##                                                               ##
##       Human cell landscape data pre-processing                ##
##                                                               ##
###################################################################

## Loading data
hcl <- readRDS('data/HCL_cell_annotations_all_seu.rds')

## Log2 normalization
hcl <- NormalizeData(hcl)
hcl <- ScaleData(hcl)

## Selecting Adult tissue
hcl$'adult' <- ! grepl('fetal', tolower(hcl$celltype)) 

saveRDS(subset(hcl, adult == TRUE),
        'data/hcl_normalized_adult_seu.rds')

###################################################################
##                                                               ##
##               Liao M. et al, 2020 pre-processing              ##
##                                                               ##
###################################################################

## Loading patient 141 mild group data
mild.h5 <- H5Fopen('data/GSM4339769_C141_filtered_feature_bc_matrix.h5') 
str(mild.h5$'matrix'$barcodes)
mild.mtx <- sparseMatrix(i = mild.h5$'matrix'$indices,
                         p = mild.h5$'matrix'$indptr, 
                         x = mild.h5$'matrix'$data,
                         dims = mild.h5$'matrix'$shape)
rownames(mild.mtx) <- mild.h5$'matrix'$features$name
colnames(mild.mtx) <- mild.h5$'matrix'$barcodes

mild.seu <- CreateSeuratObject(counts = mild.mtx,
                               project = 'sars.cov2', 
                               assay = 'RNA', 
                               min.cells = 1, 
                               min.features = 1
)

#############################################################################
## Mild cases group 142
mild.142.h5 <- H5Fopen('data/GSM4339770_C142_filtered_feature_bc_matrix.h5') 

mild.142.mtx <- sparseMatrix(i = mild.142.h5$'matrix'$indices,
                             p = mild.142.h5$'matrix'$indptr, 
                             x = mild.142.h5$'matrix'$data,
                             dims = mild.142.h5$'matrix'$shape)
rownames(mild.142.mtx) <- mild.142.h5$'matrix'$features$name
colnames(mild.142.mtx) <- mild.142.h5$'matrix'$barcodes

mild.142.seu <- CreateSeuratObject(counts = mild.142.mtx,
                                   project = 'sars.cov2', 
                                   assay = 'RNA', 
                                   min.cells = 1, 
                                   min.features = 1
)

############################################################################

## Severe group 143
sev.143.h5 <- H5Fopen('data/GSM4339771_C143_filtered_feature_bc_matrix.h5') 

sev.143.mtx <- sparseMatrix(i = sev.143.h5$'matrix'$indices,
                            p = sev.143.h5$'matrix'$indptr, 
                            x = sev.143.h5$'matrix'$data,
                            dims = sev.143.h5$'matrix'$shape)
rownames(sev.143.mtx) <- sev.143.h5$'matrix'$features$name
colnames(sev.143.mtx) <- sev.143.h5$'matrix'$barcodes

sev.143.seu <- CreateSeuratObject(counts = sev.143.mtx,
                                  project = 'sars.cov2', 
                                  assay = 'RNA', 
                                  min.cells = 1, 
                                  min.features = 1
)

#############################################################################
## Loading mild group 144
mild.144.h5 <- H5Fopen('data/GSM4339772_C144_filtered_feature_bc_matrix.h5') 

mild.144.mtx <- sparseMatrix(i = mild.144.h5$'matrix'$indices,
                             p = mild.144.h5$'matrix'$indptr, 
                             x = mild.144.h5$'matrix'$data,
                             dims = mild.144.h5$'matrix'$shape)
rownames(mild.144.mtx) <- mild.144.h5$'matrix'$features$name
colnames(mild.144.mtx) <- mild.144.h5$'matrix'$barcodes

mild.144.seu <- CreateSeuratObject(counts = mild.144.mtx,
                                   project = 'sars.cov2', 
                                   assay = 'RNA', 
                                   min.cells = 1, 
                                   min.features = 1
)

############################################################################
## Severe group 145

## Loading 145 severe case
sev.145.h5 <- H5Fopen('data/GSM4339773_C145_filtered_feature_bc_matrix.h5') 

sev.145.mtx <- sparseMatrix(i = sev.145.h5$'matrix'$indices,
                            p = sev.145.h5$'matrix'$indptr, 
                            x = sev.145.h5$'matrix'$data,
                            dims = sev.145.h5$'matrix'$shape)
rownames(sev.145.mtx) <- sev.145.h5$'matrix'$features$name
colnames(sev.145.mtx) <- sev.145.h5$'matrix'$barcodes

sev.145.seu <- CreateSeuratObject(counts = sev.145.mtx,
                                  project = 'sars.cov2', 
                                  assay = 'RNA', 
                                  min.cells = 1, 
                                  min.features = 1
)

#######################################################################


############################################################################
## Loading 146 group cases
severe.h5 <- H5Fopen('data/GSM4339774_C146_filtered_feature_bc_matrix.h5') 

severe.mtx <- sparseMatrix(i = severe.h5$'matrix'$indices,
                           p = severe.h5$'matrix'$indptr, 
                           x = severe.h5$'matrix'$data,
                           dims = severe.h5$'matrix'$shape)
rownames(severe.mtx) <- severe.h5$'matrix'$features$name
colnames(severe.mtx) <- severe.h5$'matrix'$barcodes

severe.seu <- CreateSeuratObject(counts = severe.mtx,
                                 project = 'sars.cov2', 
                                 assay = 'RNA', 
                                 min.cells = 1, 
                                 min.features = 1
)

###############################################################################
## Integration
mild.seu$orig.ident <- 'mild.141'
mild.142.seu$orig.ident <- 'mild.142'
mild.144.seu$orig.ident <- 'mild.144'
sev.145.seu$orig.ident <- 'severe.145'
severe.seu$orig.ident <- 'severe.146'
merged <- merge(x = mild.seu, y = list(mild.142.seu, mild.144.seu,
                                       sev.145.seu, 
                                       severe.seu))
rm(mild.seu, mild.142.seu, mild.144.seu, sev.145.seu, severe.seu,
   mild.h5, mild.142.h5, mild.144.h5, sev.145.h5, severe.mtx,
   mild.mtx, mild.142.mtx, mild.144.mtx, sev.145.mtx, severe.mtx)

merged <- SCTransform(merged)
merged <- FindVariableFeatures(merged, nfeatures = 1000)
merged <- RunPCA(merged)
ElbowPlot(merged)
merged <- RunUMAP(merged, dims = 1:10)
merged <- FindNeighbors(merged, dims = 1:10)
merged <- FindClusters(merged, resolution = 1.2)

DimPlot(merged, group.by = 'orig.ident', label = TRUE) + NoLegend()
DimPlot(merged, group.by = 'seurat_clusters', label = TRUE) + NoLegend()
merged$'patient_type' <- ifelse(grepl('mild', merged$orig.ident),
                                'Mild', 'Severe') 
merged$tag <- paste(merged$patient_type, merged$seurat_clusters, sep = '.')
tw_DotPlot(merged, 
           features = c('KRT18', 'TPPP3', 'SPINT2'), 
           group_by = 'seurat_clusters',
           dot_scale = 12)
FeaturePlot(merged, features = 'KRT18')

merged %>% 
        subset(seurat_clusters == 18) %>%
        tw_DotPlot(features = markers, 
                   group_by = 'patient_type',
                   cols = c('gray', 'red'), 
                   dot_scale = 20)

saveRDS(merged, '~/sc/sars-cov2/data/liao2020/liao2020_merged_seu.rds')
saveRDS(subset(merged, seurat_clusters == 18), 
        '~/sc/sars-cov2/data/liao2020/epithelial_spint2_cells.rds')

#################################################################
