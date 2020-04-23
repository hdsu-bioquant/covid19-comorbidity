## Laughney data analysis

## dependencies
library(Seurat)
library(dplyr)
library(randomForest)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(miceadds)

## Random forest for feature selection 
rforestImportance <- function(seurat, features){
        markers.vals <- FetchData(seurat, vars = features)
        
        output <- seurat@meta.data$normal_tumor
        
        rforest <- randomForest(x = markers.vals, 
                                y = as.factor(output), 
                                mtry = 5, 
                                ntree = 100)
        
        imp.df <- rforest$importance %>% 
                as.data.frame() %>% 
                mutate(gene=rownames(rforest$importance)) %>% 
                arrange(desc(MeanDecreaseGini)) 
        
        return(imp.df)
} 

## subsetting seurat inside an apply loop
subset_cell_type <- function(seurat, cell_type) {
        set <- seurat$cell_type == cell_type & seurat$normal_tumor != "MET"
        set <- (1:ncol(seurat))[set]
        seurat_by_cell_type <- subset(seurat, cells = set)
        return(seurat_by_cell_type)
}

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

##############################################################################################
## Reading Laughney dataset in seurat format
path2data <- '/media/ag-cherrmann/cramirez/covid19-comorbidity/data/'

laughney_seurat_path <- paste0(path2data, 'laughney2020/laughney_h5_norm_seurat.rds')
laughney_seu <- readRDS(laughney_seurat_path)
laughney_seu$'normal_vs_tumor' <- plyr::mapvalues(laughney_seu$normal_tumor,
                                                  from = unique(laughney_seu$normal_tumor),
                                                  to = c('METASTASIS','NORMAL', 'TUMOR'))

##############################################################################################
## Possible coreceptor of SARS-COV2

## interactome
interactome <- read.table(paste0(path2data, 'interactome.tsv'),
                          header = TRUE, 
                          stringsAsFactors = FALSE)
interactome <- interactome$PreyGene

## transcriptome
transcriptome <- read.table(paste0(path2data, 'blanco_melo_transcriptome.DEG.tsv'),
                            header = TRUE, 
                            stringsAsFactors = FALSE)
transcriptome <- transcriptome$GeneName

## Biogrid
biogrid <- read.table(paste0(path2data, 'BIOGRID-CORONAVIRUS-3.5.184.tab3.txt'),
                      header = TRUE, 
                      stringsAsFactors = FALSE, 
                      sep = '\t')
biogrid <- biogrid$Official.Symbol.Interactor.B

## Ashwin markers
ashwin <- c('ANPEP', 'DPP4', 'IFITM3', 
            'B2M', 'RPL36', 'PABPC1', 'RHOA')

## Union of all markers
markers <- unique(c(interactome, transcriptome, biogrid, ashwin,'ACE2', 'TMPRSS2'))


#########################################################################################
## Figure 1. Expression of SARS-COV2 related markers in normal and tumor lung tissues
## by cell type

## Subsetting Laughney dataset
selected_cell_types <- c('ENDOTHELIAL', 'EPITHELIAL', 
                         'FIBROBLAST', 'PERICYTE',
                         "PROLIFERATING MESENCHYMAL PROGENITOR") 

## Performing random forest for feature selection by cell type
seurat_list <- lapply(selected_cell_types, 
                  function(type) subset_cell_type(laughney_seu, type)) 
names(seurat_list) <- selected_cell_types
forests <- lapply(seurat_list, function(seurat) rforestImportance(seurat, markers))
names(forests) <- selected_cell_types

## Extracting expression values for SARS-COV2 possible coreceptors
for (i in 1:length(seurat_list)){
        mtx <- FetchData(seurat_list[i][[1]], 
                         vars = forests[i][[1]]$gene[1:50]) %>%
                as.matrix() %>% t
        #metadata <- laughney_subset@meta.data
        
        ## Heatmap of gene expression by cell type
        ## assigning colors to columns
        normal_vs_tumor_cols <- c('NORMAL' = 'green', 'TUMOR' = 'red')
        col_anns <- HeatmapAnnotation(
                Sample_type = seurat_list[i][[1]]@meta.data$normal_vs_tumor,
                col = list(Sample_type = normal_vs_tumor_cols)
        )
        H <- Heatmap(mtx, 
                top_annotation = col_anns, 
                show_column_names = FALSE,
                column_title = names(seurat_list)[i]
        )
        print(H)
}

##############################################################################
## Figure 2. Expression of top 50 SARS-COV2 associated markers using random
## forest in lung tissue

## selecting cell subtypes
laughney_subset <- laughney_seu %>% 
        subset(cell_type %in% selected_cell_types & 
                       normal_vs_tumor  != 'METASTASIS') 
laughney_subset$'abbv' <- plyr::mapvalues(laughney_subset$cell_type,
                                          from = unique(laughney_subset$cell_type),
                                          to = c("EPITHELIAL", "FIBROBLAST",
                                                 "PMP",
                                                 "ENDOTHELIAL", "PERICYTE"))

DimPlot(laughney_subset, group.by = 'abbv', label = T) + NoLegend()

## Performing random forest
forest <- rforestImportance(laughney_subset, markers)

## Extracting top 50 feature selected vars
mtx <- FetchData(laughney_subset, 
                 vars = forest$gene[1:50]) %>%
        as.matrix() %>% t
metadata <- laughney_subset@meta.data

## Plotting heatmap
normal_vs_tumor_cols <- c('NORMAL' = 'green', 'TUMOR' = 'red')
cell_type_cols <- c("ENDOTHELIAL" = 'salmon', "EPITHELIAL" = 'steelblue',
               "FIBROBLAST" = 'purple', "PERICYTE" = 'aquamarine3',
               "PMP" = 'chartreuse3')
col_anns <- HeatmapAnnotation(
        Sample_type = laughney_subset@meta.data$normal_vs_tumor,
        cell_type = laughney_subset@meta.data$abbv,
        col = list(Sample_type = normal_vs_tumor_cols, 
                   cell_type = cell_type_cols))
H <- Heatmap(mtx, 
             top_annotation = col_anns, 
             show_column_names = FALSE,
             column_title = 'Lung Tissue'
)
H

#############################################################################
## Figure 3. Dot plot showing gene expression by tissue and sample type
laughney_subset$'cell_type_normal_tumor' <- with(laughney_subset@meta.data, 
                                                 paste(abbv, 
                                                       normal_vs_tumor))
p1 <- tw_DotPlot(laughney_subset, 
        group_by = 'cell_type_normal_tumor',
        features = forest$gene[1:30], dot_scale = 8)
## change order of categories
p1$data$id <- factor(p1$data$id, levels = c("PMP NORMAL", 
                                            "ENDOTHELIAL NORMAL", 
                                            "EPITHELIAL NORMAL", 
                                            "FIBROBLAST NORMAL",
                                            "PERICYTE NORMAL", 
                                            "PMP TUMOR",
                                            "ENDOTHELIAL TUMOR", 
                                            "EPITHELIAL TUMOR",
                                            "FIBROBLAST TUMOR", 
                                            "PERICYTE TUMOR"))
plot(p1)

###############################################################################
## Figure 4. Gene expression is not biased by differences in read counts

## beta actin expression
VlnPlot(laughney_subset, 
        group.by = 'abbv',
        split.by = 'normal_vs_tumor',
        features = 'ACTB')
## umap plot showing beta actin expression
FeaturePlot(laughney_subset, 
            features = 'ACTB', 
            split.by = 'normal_vs_tumor')

################################################################################
## Fig. 5. DEG analysis normal vs tumor tissues

Idents(laughney_subset) <- laughney_subset$normal_vs_tumor
lung.markers <- FindAllMarkers(laughney_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
lung.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

## Extracting top 50 feature selected vars
deg_mtx <- FetchData(laughney_subset, 
                 vars = lung.markers$gene) %>%
        as.matrix() %>% t
metadata <- laughney_subset@meta.data

## Plotting heatmap
normal_vs_tumor_cols <- c('NORMAL' = 'green', 'TUMOR' = 'red')
cell_type_cols <- c("ENDOTHELIAL" = 'salmon', "EPITHELIAL" = 'steelblue',
                    "FIBROBLAST" = 'purple', "PERICYTE" = 'aquamarine3',
                    "PMP" = 'chartreuse3')
col_anns <- HeatmapAnnotation(
        Sample_type = laughney_subset@meta.data$normal_vs_tumor,
        cell_type = laughney_subset@meta.data$abbv,
        col = list(Sample_type = normal_vs_tumor_cols, 
                   cell_type = cell_type_cols))
H <- Heatmap(deg_mtx, 
             top_annotation = col_anns, 
             show_column_names = FALSE,
             column_title = 'Lung Tissue - DEG', 
             show_row_names = FALSE
)
H

##########################################################################
## Fig 6. Intersection of DEG and possible co-receptors for SARS-COV2
## in lung tissue

deg_int_markers <- intersect(markers, lung.markers$gene)

## Extracting top 50 feature selected vars
deg_int_markers_mtx <- FetchData(laughney_subset, 
                     vars = deg_int_markers) %>%
        as.matrix() %>% t
metadata <- laughney_subset@meta.data

gene.metadata <- data.frame(
        gene = rownames(laughney_subset),
        interactome = sapply(rownames(laughney_subset), function(g)
                ifelse( g %in%  interactome, 1, 0)),
        transcriptome = sapply(rownames(laughney_subset), function(g)
                ifelse( g %in%  transcriptome, 1, 0)),
        biogrid = sapply(rownames(laughney_subset), function(g)
                ifelse( g %in%  biogrid, 1, 0))
)
H <- Heatmap(deg_int_markers_mtx, 
             top_annotation = col_anns, 
             show_column_names = FALSE,
             column_title = 'Lung Tissue - intersect(SARS-CVO2 markers, DEG)'
)
H

##########################################################################
## Figure 7. Intersection of DEG and possible co-receptors for SARS-COV2
## by cell type

for (i in 1:length(seurat_list)) { Idents(seurat_list[i][[1]]) <- seurat_list[i][[1]]$normal_vs_tumor}
cell_type.markers <- lapply(seurat_list, function(x) 
                               FindAllMarkers(x, 
                                              only.pos = TRUE, 
                                              min.pct = 0.25, 
                                              logfc.threshold = 0.25))
cell_type.markers <- lapply(cell_type.markers, function(x)
                                 x %>% group_by(cluster) %>% 
                                    top_n(n = 10, wt = avg_logFC))

names(cell_type.markers) <- selected_cell_types
for (i in 1:length(cell_type.markers)) {
        cell_type.markers[i][[1]]$'cell_type' <- names(cell_type.markers[i])
}

lapply(cell_type.markers, function(x) intersect(x$gene, markers))

tw_DotPlot(laughney_subset, 
        features = c("CXCL2", "CXCL1", ## Downregulated in epithelial tumor sample
                     "IFI27",          ## Downregulated in fibroblast tumor sample
                     'ACTB',           ## reference house keeping gene
                     'ACE2', 'TMPRSS2'),   ## reference SARS-COV2 markers 
        group_by = 'cell_type_normal_tumor')


cell_type.markers.merged <- bind_rows(cell_type.markers)
cell_type.markers.merged$'int' = sapply(cell_type.markers.merged$gene, 
                                                      function(x) 
                                                              ifelse(x %in% markers,
                                                                     17, 19))
cell_type.markers.merged$'gene_lab' = sapply(cell_type.markers.merged$gene, 
                                        function(x) 
                                                ifelse(x %in% markers,
                                                       x, ''))
ggplot(cell_type.markers.merged, aes(avg_logFC, -log10(p_val_adj),
                                        colour = cell_type)) + 
                                geom_point(aes(shape = int)) +
               geom_text(aes(label=gene_lab),hjust=0, vjust=0) +
        scale_shape_identity() + theme_classic()

path2analysis <- '/media/ag-cherrmann/cramirez/covid19-comorbidity/'
write.table(cell_type.markers.merged, 
            paste0(path2analysis, 'top10_DEG_lung_by_cell_type.tsv'),
            row.names = FALSE, 
            sep = '\t')

##########################################################################
## CXCL1 and -2 are expressed in normal lung tissue in an independet dataset
hcl <- readRDS('~/Downloads/HCL_lung_seurat.rds')
hcl <- subset(hcl, celltype %in% c('Endothelial cell (endothelial to mesenchymal transition)',
                                   'Endothelial cell (APC)',
                                   'Epithelial cell (intermediated)',
                                   'Endothelial cell',
                                   'Epithelial cell',
                                   'Fibroblast'))
hcl <- NormalizeData(hcl) %>% ScaleData()
hcl <- FindVariableFeatures(hcl, nfeatures = 1000)
hcl <- RunPCA(hcl)
hcl <- RunUMAP(hcl, dims = 1:20)
hcl <- NormalizeData(hcl) %>% ScaleData()
DotPlot(hcl, features = c('CXCL1', 'CXCL2', 'IFI27'), group.by = 'celltype')
DimPlot(hcl, group.by = 'celltype', label = TRUE) + NoLegend()
FeaturePlot(hcl, features = c('CXCL1', 'CXCL2', 'IFI27'))
