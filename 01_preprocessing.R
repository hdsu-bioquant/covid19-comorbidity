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

##################################################################
##                                                              ##
##              Processing Colorectal Cancer scRNA-Seq          ##
##      Reference: https://www.nature.com/articles/ng.3818      ##
##                                                              ##
##################################################################
## Downloading data
url1 <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81861&format=file&file=GSE81861%5FCRC%5FNM%5Fall%5Fcells%5FCOUNT%2Ecsv%2Egz'
url2 <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81861&format=file&file=GSE81861%5FCRC%5FNM%5Fall%5Fcells%5FFPKM%2Ecsv%2Egz'
url3 <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81861&format=file&file=GSE81861%5FCRC%5FNM%5Fepithelial%5Fcells%5FCOUNT%2Ecsv%2Egz'
url4 <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81861&format=file&file=GSE81861%5FCRC%5FNM%5Fepithelial%5Fcells%5FFPKM%2Ecsv%2Egz'
url5 <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81861&format=file&file=GSE81861%5FCRC%5Ftumor%5Fall%5Fcells%5FCOUNT%2Ecsv%2Egz'
url6 <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81861&format=file&file=GSE81861%5FCRC%5Ftumor%5Fall%5Fcells%5FFPKM%2Ecsv%2Egz'
url7 <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81861&format=file&file=GSE81861%5FCRC%5Ftumor%5Fepithelial%5Fcells%5FCOUNT%2Ecsv%2Egz'
url8 <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81861&format=file&file=GSE81861%5FCRC%5Ftumor%5Fepithelial%5Fcells%5FFPKM%2Ecsv%2Egz'
url9 <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81861&format=file&file=GSE81861%5FCell%5FLine%5FCOUNT%2Ecsv%2Egz'
url10 <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81861&format=file&file=GSE81861%5FCell%5FLine%5FFPKM%2Ecsv%2Egz'
url11 <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81861&format=file&file=GSE81861%5FGEO%5FEGA%5FID%5Fmatch%2Ecsv%2Egz'

urls <- c(url1, url2, url3, url4, url5, url6, url7, url8, url9, url10, url11)

## Downloading files
for (u in urls) {
        print(u)
        file.name <- gsub('.*file=', '', u)
        file.name <- gsub('%', '_', file.name)
        file.name <- gsub('_2Ecsv_2Egz', '.csv.gz', file.name)
        print(file.name)
        download.file(url = u, destfile = paste0('data/', file.name))
}

gzips <- list.files()
gzips <- gzips[!grepl('match', gzips)]
gzips

sapply(gzips, function(g) gunzip(g, remove=FALSE))

#################################################################################################################
## Processing normal epithelial cells
normal_epi_counts <- read.csv("GSE81861_5FCRC_5FNM_5Fepithelial_5Fcells_5FCOUNT.csv")
dim(normal_epi_counts)

## Processing gene names
gene <- sapply(normal_epi_counts$X, function(x) { 
        sp_name <- str_split(x, pattern = '_')[[1]][2]
}) 
normal_epi_counts <- select(normal_epi_counts, -X)
## removing duplicated ges
dups <- duplicated(gene)
normal_epi_counts <- normal_epi_counts[!dups,]
rownames(normal_epi_counts) <- gene[!dups] 


## Extracting annotations from column
## cell type
cell_type <- sapply(colnames(normal_epi_counts), function(x) { 
        sp_name <- str_split(x, pattern = '__')[[1]][2]
}) 
sample_ann1 <- sapply(colnames(normal_epi_counts), function(x) { 
        sp_name <- str_split(x, pattern = '__')[[1]][1]
}) 
sample_ann2 <- sapply(colnames(normal_epi_counts), function(x) { 
        sp_name <- str_split(x, pattern = '__')[[1]][3]
})
anns <- data.frame('cell_type' = cell_type,
                   'sample_ann1' = sample_ann1,
                   'sample_ann2' = sample_ann2)

## Definition of the seurat object
epi_normal_seu <- CreateSeuratObject(
        counts = normal_epi_counts,
        project = 'sars-cov2',
        assay = 'RNA',
        min.cells = 1,
        min.features = 1
)
epi_normal_seu <- AddMetaData(epi_normal_seu, 
                              metadata = anns)
epi_normal_seu$'normal_vs_tumor' <- 'normal'

################################################################################################################
## Processing tumor epithelial cells

tumor_epi_counts <- read.csv("GSE81861_5FCRC_5Ftumor_5Fepithelial_5Fcells_5FCOUNT.csv")
dim(tumor_epi_counts)
tumor_epi_counts[1:5, 1:5]

## Processing gene names
gene <- sapply(tumor_epi_counts$X, function(x) { 
        sp_name <- str_split(x, pattern = '_')[[1]][2]
})

## Removing gene names from the count matrix
tumor_epi_counts <- select(tumor_epi_counts, -X)

## removing duplicated ges
dups <- duplicated(gene)
tumor_epi_counts <- tumor_epi_counts[!dups,]
rownames(tumor_epi_counts) <- gene[!dups] 


## Extracting annotations from column
## cell type
cell_type <- sapply(colnames(tumor_epi_counts), function(x) { 
        sp_name <- str_split(x, pattern = '__')[[1]][2]
}) 
sample_ann1 <- sapply(colnames(tumor_epi_counts), function(x) { 
        sp_name <- str_split(x, pattern = '__')[[1]][1]
}) 
sample_ann2 <- sapply(colnames(tumor_epi_counts), function(x) { 
        sp_name <- str_split(x, pattern = '__')[[1]][3]
})
anns <- data.frame('cell_type' = cell_type,
                   'sample_ann1' = sample_ann1,
                   'sample_ann2' = sample_ann2)

## Definition of the seurat object
epi_tumor_seu <- CreateSeuratObject(
        counts = tumor_epi_counts,
        project = 'sars-cov2',
        assay = 'RNA',
        min.cells = 1,
        min.features = 1
)
epi_tumor_seu <- AddMetaData(epi_tumor_seu, 
                             metadata = anns)
epi_tumor_seu$'normal_vs_tumor' <- 'tumor'

#############################################################################################################
counts_vln <- VlnPlot(epi, features = 'nCount_RNA', group.by = 'normal_vs_tumor')
counts_umap <- FeaturePlot(epi, features = 'nCount_RNA', pt.size = 3)

## Mitochondrial percentage
epi[["percent.mt"]] <- PercentageFeatureSet(epi, pattern = "^MT-")
mito_vln <- VlnPlot(epi, features = 'percent.mt', group.by = 'normal_vs_tumor')
mito_umap <- FeaturePlot(epi, features = 'percent.mt', pt.size = 3)

## Quality plots
CombinePlots(list(counts_vln, mito_vln, counts_umap, mito_umap))

#############################################################################################################
## Merge, normalization and explorative analysis
epi <- merge(x = epi_normal_seu, y = epi_tumor_seu)

epi <- SCTransform(epi)
epi <- FindVariableFeatures(epi, nfeatures = 3000)
epi <- RunPCA(epi, features = VariableFeatures(epi))
ElbowPlot(epi)
epi <- RunUMAP(epi, dims = 1:20)
epi <- FindNeighbors(epi)
epi <- FindClusters(epi, resolution = 1)

## Saving integrated data
saveRDS(epi, 'data/GSE81861_CRC/epi_seu.rds')

###################################################################
##                                                               ##
##              Behjati et al, 2018 data                         ##
##              renal Clear Cell Carcinoma data                  ##
##                                                               ##
##################################################################


## Count matrix
kidney <- readMM('data/tableOfCounts.mtx')
## gene annotations
gene_anns <- read.table('data/aat1699_DataS1/tableOfCounts_rowLabels.tsv',
                        header = T)
## cell annotations
anns <- readRDS('data/cell_annot.RDS')
barcodes <- read.table('data/tableOfCounts_colLabels.tsv',
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
saveRDS(behjati, 'data/behjati2018_seu.rds')

####################################################################
##                                                                ##
##              Hepatocyte Cell Carcinoma scRNA-Seq data          ##
##                                                                ##
####################################################################

markers <- c('KRT18', 'KRT7', 'KRT19', 'EPCAM', 'SOX9', ## Endothelial cells
             'ALB', 'TF', 'CYP3A4', 'CYP2E1', 'APOB', 'ASGR1', 'PCK1', 'APOE',  ## Hepatocytes
             'CLEC4G', 'CLEC4M', 'FLT1', 'PECAM1',  ## Liver Sinusoidal
             'KLRB1', 'PTPRC', 'CD3E',    ## Immune cells
             'CD163', 'MAFB', 'VSIG4'    ## Kupffer cells
)

## Lu data on Hepatocellular carcinoma scRNA-Seq
url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE149614&format=file&file=GSE149614%5FHCC%2EscRNAseq%2ES71915%2Ecount%2Etxt%2Egz'

download.file(url = url, destfile = 'data/GSE149614_HCC.scRNAseq.S71915.count.txt.gz')

## Reading data
liver.mtx <- read.table('data/GSE149614_HCC.scRNAseq.S71915.count.txt.gz', 
                        sep = '\t', quote = '', header = TRUE)
dim(liver.mtx)

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

## Clustering
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

saveRDS(liver_seu.ss, 'data/GSE149614_liver_seu.rds')

