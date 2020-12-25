## Evaluation of the permissivity signature in mild vs severe
## cases

## Dependencies
library(ggplot2)
library(dplyr)
library(ggrepel)
library(reshape)
library(maditr)
library(pheatmap)
library(biomaRt)
library(Seurat)
library(RColorBrewer)
library(cmapR)
library(ggpubr)
library(viridis)
library(ggradar)


set.seed(333)

## Setting up paths
project_dir <- '/Users/carlosramirez/sc/covid19-comorbidity/'
setwd(project_dir)

#####################################################################
##                                                                 ##
##      Figure 1a. TF bound to SPINT2 and TMPRSS2 genomic          ##
##      loci                                                       ##
##                                                                 ##
#####################################################################

## Loading TOBIAS results
files <- list.files('analysis/TOBIAS/tfbs_HIO/', 
                    full.names = TRUE)
tfbs.list <- lapply(files,
                    function(x)
                            read.table(x, header = TRUE))

## Merging TFBS
tfbs.df  <- do.call(rbind, tfbs.list)

## Extracting binding sites and restricting to SPINT2
## and TMPRSS2 sites
binding_sites <- filter(tfbs.df, 
                        gene_name %in% c('SPINT2', 'TMPRSS2',
                                         'SPINT1', 'ST14'))
binding_sites <- mutate(binding_sites, 
                        TFBS_name= gsub('_.*', '', TFBS_name))

## Restricting to positively correlated TFs in the TF activities
## scRNA-Seq analysis
pos.cor <- read.table('~/sc/covid19-comorbidity/analysis/network_tobias.tsv',
                      header = TRUE, 
                      sep = '\t')
tfbs <- pos.cor$TFBS_name %>% unique()

binding_sites <- filter(binding_sites, 
                        TFBS_name %in% tfbs)
head(binding_sites)
dim(binding_sites)

## Creating list of gene coordinates
ensembl <- useMart("ensembl", 
                   dataset="hsapiens_gene_ensembl")
values <- c('SPINT2', 'TMPRSS2',
            'SPINT1', 'ST14')
gene.coor <- getBM(attributes= c('chromosome_name',
                                 'start_position',
                                 'end_position',
                                 'hgnc_symbol'), 
                   filters = 'hgnc_symbol', 
                   values = values, 
                   mart = ensembl)
gene.coor <- mutate(gene.coor, 
                    chromosome_name = paste0('chr',
                                             chromosome_name))
gene.coor <- mutate(gene.coor,
                    start_position = start_position - 5000,
                    end_position = end_position + 5000)

## Saving the table
write.table(gene.coor,
            file = 'analysis/gene_coordinates.tsv',
            sep = '\t',
            col.names = FALSE,
            row.names = FALSE, 
            quote = FALSE)

#######################################################################
##                                                                   ##
##      Figure 1B. Transcription factor activities correlations      ##
##      to SPINT2 or TMPRSS2 gene expression                         ##
##                                                                   ##
#######################################################################

## Loading organoids scRNA-Seq dataset
load('~/sc/organoids/data/COVID19_July.rda')

## Loading TF activities
TF_act <- read.csv(
        '~/sc/organoids/analysis/scenic/ileum_corrected/ileum_aucell.csv'
)

## Calculation of the mean expression values
tfa_mean <- apply(dplyr::select(TF_act, -Cell), 
                  2, mean)
names(tfa_mean) <- gsub('\\.', '', names(tfa_mean))

## mean gene expression
tfs <- gsub('\\.', '', colnames(dplyr::select(TF_act, -Cell)))
tfs.mean_exp <- Illeum_H_T %>%
        FetchData(vars = tfs) %>%
        apply(2, mean)

## tfa correlation vs gene expression in ileum
spint2_cor_tfa <- ileum.cor['SPINT2', ]
names(spint2_cor_tfa) <- gsub('\\.', '', names(spint2_cor_tfa)) 
any(names(spint2_cor_tfa) != names(tfs.mean_exp) )
ileum.tfa.df <- data.frame(
        spint2_cor_tfa = spint2_cor_tfa,
        mean_expression = tfs.mean_exp
)
dim(ileum.tfa.df)

## Adding TMPRSS2 coregulation
tmprss2_cor_tfa <- ileum.cor['TMPRSS2', ]
ileum.tfa.df$'tmprss2_cor_tfa' <- tmprss2_cor_tfa

st14_cor_tfa <- ileum.cor['ST14', ]
ileum.tfa.df$'st14_cor_tfa' <- st14_cor_tfa

spint1_cor_tfa <- ileum.cor['SPINT1', ]
ileum.tfa.df$'spint1_cor_tfa' <- spint1_cor_tfa

## Adding mean TF activities
ileum.tfa.df$'TF_activity_mean' <- 0
ileum.tfa.df[names(tfa_mean),]$'TF_activity_mean'  <- tfa_mean

markers <- c('ONECUT3',
             'IRF1',
             'IRF7',
             'FOXC1',
             'JUNB',
             'JUND',
             'FOS', 
             'FOSL1',
             'ELF3',
             'KLF4')

pdf('figures/tfa_cor_spint2_tmprss2_ileum_org_v2.pdf',
    height = 7, width = 8)
ileum.tfa.df %>%
        mutate(gene = rownames(ileum.tfa.df)) %>%
        mutate(highlight = ifelse( gene %in% markers,
                                   TRUE, FALSE)) %>%
        mutate(gene_label = ifelse(highlight == TRUE,
                                   gene, '')) %>%
        ggplot(aes(x=spint2_cor_tfa, 
                   y=tmprss2_cor_tfa,
                   colour=mean_expression,
                   label = gene_label,
                   size = TF_activity_mean)) +
        geom_point() +
        geom_text_repel(colour='black',
                        size=4) +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        geom_hline(yintercept = cor_th, 
                   linetype='dashed') +
        geom_vline(xintercept = cor_th,
                   linetype = 'dashed') +
        scale_color_viridis() +
        ylab('Cor(TF activity, TMPRSS2)') +
        xlab('Cor(TF activity, SPINT2)') +
        ggtitle('Ileum organoid')
dev.off()

#######################################################################
##                                                                   ##
##      Figure 1c. Network of regulatory interactions from TFs       ##
##      to SPINT2 and TMPRSS2 as inferred by TF footprinting         ##
##      analysis                                                     ##
##                                                                   ##
#######################################################################

## Saving results to tsv file for importing into cytoscape
write.table(binding_sites,
            file = 'analysis/network_spint2_tmprss2.tsv',
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t',
            quote = FALSE)

#######################################################################
##                                                                   ##
##      Figure 1D. SPINT2/TMPRSS2 correlations across tissues        ##
##      using GTEx data                                              ##
##                                                                   ##
#######################################################################

url_file <- 'https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz'
file_path <- 'data/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct'

## Downloading 
if ( ! file.exists(file_path) ){
        download.file(url = url_file,
                      destfile = file_path)
}

gtex <- parse_gctx(file_path)
dim(gtex@mat)
gtex@mat[1:5, 1:5]
saveRDS(gtex@mat, 'data/gtex/GTEx_tpm.mtx')
rm(gtex)

gtex.mtx <- readRDS('data/gtex/GTEx_tpm.mtx')

## Gene annotation
## Reading gene annotations
gtf_url <- 'https://storage.googleapis.com/gtex_analysis_v8/reference/gencode.v26.GRCh38.genes.gtf'
#download.file(gtf_url, destfile = 'data/gtex/gencode.v26.GRCh38.genes.gtf')

gtf <- read.table(file = gtf_url, sep = '\t')
head(gtf)


extract_pattern <- function(string_vec, pattern){
        string_vec <- gsub(paste0('.*', pattern), '', string_vec)
        string_vec <- gsub(' |;.*', '', string_vec)
        
        return(string_vec)
}

## Processing GTEx gene annotations 
gene_id <- extract_pattern(gtf$V9, 'gene_id')
gene_name <- extract_pattern(gtf$V9, 'gene_name')

gene_anns <- data.frame(gene_id=gene_id,
                        gene_name=gene_name)
gene_anns <- gene_anns[!duplicated(gene_anns), ]
rownames(gene_anns) <- gene_anns$gene_id

## Cleanning gene names
gtex.mtx <- gtex.mtx[rownames(gtex.mtx) %in% 
                             rownames(gene_anns), ] 
gtex.mtx <- gtex.mtx[rownames(gene_anns), ]
rownames(gtex.mtx) <- gene_anns$gene_name
gtex.mtx[1:5, 1:5]

saveRDS(gtex.mtx, 'data/gtex/GTEx_tpm.mtx')

#gtex.mtx <- readRDS('data/gtex/GTEx_tpm.mtx')

## Sample annotations
gtex.ann <- read.table(
        'https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt',
        header = TRUE, sep = '\t',
        quote = ''
)
dim(gtex.ann)

## Subsetting to serine proteases related genes and TFs

## SPRGs
uniprot_SPRG.filtered <- readRDS('data/serine_proteases_filtered.rds')
SPRGs <- strsplit(uniprot_SPRG.filtered$Gene.names, split = ' |/') %>%
        unlist() %>% sort() %>% unique()

## Downloading list of TFs
tfs_url <- 'https://raw.githubusercontent.com/aertslab/pySCENIC/master/resources/hs_hgnc_curated_tfs.txt'
tfs <- readLines(tfs_url)

## selected markers
markers <- c('ACE2', 'CDX1', 'CDX2', 'HNF4A')

genes <- c(SPRGs, tfs, markers)

## extracting gene expression
gtex.df <- t(gtex.mtx[rownames(gtex.mtx) %in% genes, ])

## annotation by cell type
unique(gtex.ann$SMTS)
rownames(gtex.ann) <- gtex.ann$SAMPID

anns <- gtex.ann[rownames(gtex.ann) %in%
                         rownames(gtex.df), ]
## test gene expression and annotations match
any(rownames(anns) != rownames(gtex.df))

## merge ge and annotations
ge.ann <- cbind(gtex.df, anns)
ge.ann <- ge.ann[, ! duplicated(colnames(ge.ann))]

saveRDS(ge.ann, 'data/gtex/GTEx_tpm.mtx')

#ge.ann <- readRDS('data/gtex/GTEx_tpm.mtx')

selected_tissues <- c("Blood Vessel",
                      "Heart",
                      "Brain",
                      "Lung", 
                      "Pancreas",
                      "Esophagus",
                      "Stomach",
                      "Colon",
                      "Small Intestine",
                      "Prostate",
                      "Blood",
                      "Liver",
                      "Kidney")

## vector to store correlations
tissue_type <- unique(ge.ann$SMTS)
gtex.cors <- data.frame(
        tissue_type=tissue_type,
        spearman_cor=numeric(length(tissue_type)),
        p_val=numeric(length(tissue_type)),
        rsq=numeric(length(tissue_type)),
        n=numeric(length(tissue_type))
)

for (i in 1:length(tissue_type)){
        tissue <- filter(ge.ann, SMTS==tissue_type[i])
        if ( nrow(tissue) > 0 ) {
                lmod <- lm(TMPRSS2 ~ SPINT2, data = tissue)
                gtex.cors$spearman_cor[i] <- lmod$coefficients[2]
                gtex.cors$rsq[i] <- summary(lmod)$adj.r.squared
                gtex.cors$p_val[i] <- summary(lmod)$coefficients[2,4]
                gtex.cors$n[i] <- nrow(tissue)
        }
}

head(gtex.cors)

##############################################################
##                                                          ##
##      Correlation of TMPRSS2 and SPINT2 in Nasopharynx    ##
##                                                          ##
##############################################################

mtx.df <- readRDS('analysis/nasopharynx_annotated_mtx.rds')

lmod <- lm(TMPRSS2 ~ SPINT2, data = mtx.df)
nasopha <- data.frame(
        tissue_type = 'Nasopharynx',
        spearman_cor = lmod$coefficients[2],
        rsq = summary(lmod)$adj.r.squared,
        p_val = summary(lmod)$coefficients[2,4],
        n = nrow(mtx.df)
)

gtex.cors <- rbind(gtex.cors, 
                   nasopha)

gtex.cors <- filter(gtex.cors, 
                    tissue_type %in% c('Nasopharynx',
                                       selected_tissues))

rownames(gtex.cors) <- gtex.cors$tissue_type 
radar.df <- gtex.cors %>%
        arrange(desc(rsq)) %>%
        select(rsq) %>%
        t() %>%
        as.data.frame()


## Plotting
pdf('figures/spint2_tmprss2_cor_gtex_radar.pdf')
ggradar(radar.df, 
        values.radar = c(0, 0.5, 1),
        grid.min = 0, 
        plot.extent.x.sf= 1,
        plot.extent.y.sf = 1.5, 
        gridline.label.offset = -0.1, 
        group.line.width = 1,
        group.point.size = 3, 
        group.colours = 'steelblue', 
        background.circle.colour = 'white')
dev.off()

#######################################################################
##                                                                   ##
##      Figure 1E. Gene expression of ACE2, SPINT2 and TMPRSS2       ##
##      across cell types in HCL                                     ##
##                                                                   ##
#######################################################################
## Gene set visualization in Human Cell Landscape
## Reading Human Cell Landscape data
hcl <- readRDS('~/Downloads/hcl_normalized_adult_seu.rds')

## Subsetting to adult tissue
hcl$'adult' <- grepl('adult', tolower(hcl$orig.ident)) &
        ! ( grepl('fetal', tolower(hcl$celltype)))
hcl <- subset(hcl, adult == TRUE)

## Selecting cell types
celltypes <- c('Epithelial cell',
               'Endothelial cell',
               'Enterocyte progenitor',
               'Epithelial cell (intermediated)',
               'AT2 cell',
               'Goblet cell',
               'Stratified epithelial cell',
               'Kidney intercalated cell',
               'Loop of Henle',
               'Enterocyte',
               'Enterocyte progenitor',
               'Gastric chief cell',
               'Gastric endocrine cell',
               'Hepatocyte/Endodermal cell',
               'Fibroblast',
               'Macrophage',
               'M2 Macrophage',
               'Stromal cell',
               'Smooth muscle cell',
               'Endothelial cell (APC)',
               'Endothelial cell (endothelial to mesenchymal transition)',
               'Sinusoidal endothelial cell',
               'Ventricle cardiomyocyte',
               'Proximal tubule progenitor',
               'Pancreas exocrine cel',
               'Dendritic cell',
               'hESC')

markers <- c('SPINT2', 'ACE2', 'TMPRSS2')

pdf('~/sc/sars-cov2/figures/figure_01_e.pdf',
    width = 8, height = 3.5)
celltype.order <- hcl_markers$data %>% 
        subset(features.plot == 'SPINT2') %>%
        arrange(avg.exp) %>%
        select(id) %>% unlist()
celltypes <- factor(hcl_markers$data$id, levels = celltype.order)
hcl_markers$data %>%
        filter(features.plot %in% markers) %>%
        mutate(id = factor(id, levels = celltype.order[length(celltype.order):1])) %>%
        ggplot(aes(x = features.plot, 
                   y = id,
                   size = pct.exp)) +
        theme_bw(base_size = 9) +
        geom_point(aes(colour=avg.exp)) + 
        geom_point(shape = 1,colour = "black") +
        scale_size('pct.exp') +
        scale_colour_viridis() + 
        theme(axis.title = element_blank(),
              panel.grid = element_blank(), 
              axis.line = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle('Gene expression of target markers in HCL') +
        coord_flip()
dev.off()


#######################################################################
##                                                                   ##                     
##      Figure 1a. Inference of the permissivity signature (i') =    ##
##      DEG Calu vs H1299 (z) - ( DEG Calu inf vs mock (x)           ##
##                               + DEG H1299 inf vs mock (y) )       ## 
##                                                                   ##
#######################################################################


## Reading data
calu <- readRDS('data/200408.Seurat_Calu_CoV_1000_Merged.rds')
calu$strain[calu$strain == 'nan'] <- 'mock'
calu.nin <- subset(calu, orig.ident %in% c("Calu3-mock-4h-A", "Calu3-mock-4h-B") )
rm(calu)
h1299 <- readRDS('data/200406.Seurat_H1299_CoV_1000_Merged.rds')
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

## Visualizing clusters
svg(paste0(path2_sc_deg, 'calu_h1299_mock4hr_clusters.svg'))
DimPlot(calu.h1299.nin, group.by = 'orig.ident')
dev.off()

## DEG analysis
calu.h1299.nin$'cell_line' <- ifelse(grepl('Calu', calu.h1299.nin$orig.ident),
                                     'Calu', 'H1299')
Idents(calu.h1299.nin) <- calu.h1299.nin$cell_line
deg.markers <- FindAllMarkers(calu.h1299.nin)

saveRDS(deg.markers,
        paste0(path2_sc_deg,'deg_calu_h1299.rds'))

############################################################
## DEG Calu infected vs mock 12hrs
## Reading data
calu <- readRDS('~/Downloads/archive-2/200408.Seurat_Calu_CoV_1000_Merged.rds')
calu.S2.12hr.sets <- c("Calu3-mock-12h-A", "Calu3-mock-12h-B",
                       "Calu3-S2-12h-A", "Calu3-S2-12h-B")
calu.S2.12hr <- subset(calu, orig.ident %in% calu.S2.12hr.sets) 
rm(calu)

## Renaming clusters as infected and non-infected
calu.S2.12hr$'inf_vs_nin' <- ifelse(calu.S2.12hr$orig.ident %in% 
                                            c("Calu3-S2-12h-A", "Calu3-S2-12h-B"),
                                    'infected', 'non-infected') 
svg(paste0(path2_sc_deg, 'calu.inf_vs_nin.svg'))
DimPlot(calu.S2.12hr, group.by = 'inf_vs_nin', label = TRUE) + NoLegend() 
dev.off()

## Checking cell viability in Calu infected and non-Infected
mitocondrial_genes <- grep('MT-', rownames(calu.S2.12hr), value = T)
calu.S2.12hr$'mito' <- get_sign_scores(calu.S2.12hr, 
                                       gene_set = mitocondrial_genes)
VlnPlot(calu.S2.12hr, features = 'sarscov2_expression')
VlnPlot(calu.S2.12hr, features = 'mito')
VlnPlot(calu.S2.12hr, features = 'nCount_RNA')

## Filtering cell with > 50,000 counts
calu.S2.12hr <- subset(calu.S2.12hr, nCount_RNA < 50000)
VlnPlot(calu.S2.12hr, features = 'nCount_RNA')

## Normalization and dimensional reduction
calu.S2.12hr <- NormalizeData(calu.S2.12hr) %>% ScaleData()
calu.S2.12hr <- FindVariableFeatures(calu.S2.12hr, nfeatures = 3000)
calu.S2.12hr <- RunPCA(calu.S2.12hr, features = VariableFeatures(calu.S2.12hr))
ElbowPlot(calu.S2.12hr)
calu.S2.12hr <- RunUMAP(calu.S2.12hr, reduction = 'pca', dims = 1:20)

## Calu 12 hrs clusters
## **Two clearly separated clusters are observed at 12 hrs**
svg(paste0(path2_sc_deg, 'calu.mock12hr.clusters.svg'))
DimPlot(calu.S2.12hr, group.by = 'orig.ident')
dev.off()

## Checking sars-cov2 expression
## SARS-COV2 genes
scov2_genes <- grep('CoV2', rownames(calu.S2.12hr), value = T)
calu.S2.12hr$'sarscov2_expression' <- get_sign_scores(calu.S2.12hr, 
                                                      gene_set = scov2_genes)
FeaturePlot(calu.S2.12hr, features = 'sarscov2_expression')

Idents(calu.S2.12hr) <- calu.S2.12hr$inf_vs_nin
calu.signature <- FindAllMarkers(calu.S2.12hr)

################################################################################
## Getting the H1299 signature
h1299 <- readRDS('data/200406.Seurat_H1299_CoV_1000_Merged.rds')

h1299 <- NormalizeData(h1299) %>% ScaleData()
h1299 <- FindVariableFeatures(h1299, nfeatures = 1000)
h1299 <- RunPCA(h1299, features = VariableFeatures(h1299))
ElbowPlot(h1299)
h1299 <- RunUMAP(h1299, reduction = 'pca', dims=1:10)

## Visualizing conditions and clusters
svg(paste0(path2_sc_deg, 'h1299.clusters.svg'))
DimPlot(h1299, group.by = 'orig.ident')
dev.off()
svg(paste0(path2_sc_deg, 'h1299.exp.conditions.svg'))
h1299$'tag' <- sapply(h1299$strain, function(x) 
        ifelse(x == 'nan', 'Mock', x))
DimPlot(h1299, group.by = 'tag')
dev.off()

## Visualize Viral gene expression
## Viral infected cells are a combination of 12, 24 and 36 hrs treated 
h1299$'sarscov2_expression' <- get_sign_scores(h1299, gene_set = scov2_genes)
svg(paste0(path2_sc_deg, 'sarscov2_gene_expression_h1299.svg'))
FeaturePlot(h1299, features = 'sarscov2_expression')
dev.off()

## Subsetting infected cells (sars-cov2 > 0) and mock 4hrs
h1299.inf <- subset(h1299, sarscov2_expression > 0 | 
                            ( orig.ident %in% c("H1299-mock-4h-A", "H1299-mock-4h-B")))
rm(h1299)

h1299.inf <- NormalizeData(h1299.inf) %>% ScaleData()
h1299.inf <- FindVariableFeatures(h1299.inf, nfeatures = 3000)
h1299.inf <- RunPCA(h1299.inf, features = VariableFeatures(h1299.inf))
ElbowPlot(h1299.inf)
h1299.inf <- RunUMAP(h1299.inf, reduction = 'pca', dims=1:10)

svg(paste0(path2_sc_deg, 'h1299.inf.svg'))
DimPlot(h1299.inf, group.by = 'orig.ident')
dev.off()

FeaturePlot(h1299.inf, features = 'sarscov2_expression')
tw_DotPlot(h1299.inf, features = scov2_genes, 
           group_by = 'orig.ident' )

## Infected cells defined as cells with gene expression score > 0
h1299.inf$'infection_status' <- ifelse(h1299.inf$sarscov2_expression > 0,
                                       'infected', 'non-infected')

svg(paste0(path2_sc_deg, 'h1299.infection.clusters.svg'))
DimPlot(h1299.inf, group.by = 'infection_status')
dev.off()

Idents(h1299.inf) <- h1299.inf$infection_status
h1299.infection.signature <- FindAllMarkers(h1299.inf)

h1299.infection.signature %>%
        saveRDS(paste0(path2_sc_deg, 'h1299.infection.signature.rds'))

############################################################################
## Filtering induced genes upon infection
filter_degs <- filter(deg.markers, ( ! gene %in% infection.signature$gene) &
                              ( ! gene %in% h1299.infection.signature$gene) )

filter_degs %>% 
        saveRDS(paste0(path2_sc_deg, 'deg.calu_vs_h1299_mock.4hr.filtered.rds'))




##################################################################################
##                                                                              ##
##      Figure 1b. Pie chart for all the top 25% most important genes           ##
#       in various categories                                                   ##
##                                                                              ##
##################################################################################

rf <- read.table('data/variable_importance_ranked_genes.tsv',
                 header = TRUE)
head(rf)

signature.df <- read.table('data/permissivity_signature_annotated.tsv',
                           header = TRUE)
rownames(signature.df) <- signature.df$gene

signature.df <- merge(signature.df, rf)
signature.df <- arrange(signature.df, desc(IncNodePurity))
rownames(signature.df) <- signature.df$gene

## Selecting 25% top ranked genes
signature.sel <- signature.df[1:round(nrow(signature.df)*0.25), ]

signature.sel <- select(signature.sel, 
                        go_viral:infection_signature)
freqs <- apply(signature.sel, 2, sum)

names <- c("GO viral process", "Surfaceome", 
           "SARS-CoV2 Interactome", "CellphoneDB", 
           "SARS-CoV2 Infection signature")
freqs.df <- data.frame(category = names(freqs),
                       frequency = freqs,
                       abbv = names)

freqs.df <- freqs.df %>% 
        arrange(desc(category)) %>%
        mutate(prop = frequency / sum(freqs.df$frequency) * 100) %>%
        mutate(ypos = cumsum(prop)- 0.5*prop) %>%
        mutate(color = brewer.pal(nrow(freqs.df), "Set1"))

pdf('figures/pie_chart.pdf')
ggplot(freqs.df, aes(x="", y=prop, fill=abbv)) + 
        theme_void(base_size = 9) + 
        labs(fill = "") + geom_bar(stat="identity", 
                                   width=1, color="white") +
        coord_polar("y", start=0) +
        geom_text(aes(y = ypos, label = frequency), 
                  color = "white", size = 2) +
        scale_fill_manual(values = rev(freqs.df$color))
dev.off()

##################################################################
##                                                              ##
##      Figure 1b. Ranking of the genes in the permissivity     ##
##      signature                                               ##
##                                                              ##
##################################################################
calu.scov2.inf <- subset(calu, strain == 'SARSCoV2' & 
                                 infect == 'infected' &
                                 orig.ident %in% c("Calu3-S2-12h-A", 
                                                   "Calu3-S2-12h-B")
)
rm(calu)

## normalization
calu.scov2.inf <- SCTransform(calu.scov2.inf)
calu.scov2.inf <- FindVariableFeatures(calu.scov2.inf, 
                                       nfeatures = 1000)
calu.scov2.inf <- RunPCA(calu.scov2.inf, 
                         features = VariableFeatures(calu.scov2.inf)) 
ElbowPlot(calu.scov2.inf)
calu.scov2.inf <- RunUMAP(calu.scov2.inf, dims = 1:20)

## sars-cov2 expression
calu.scov2.inf$'sarscov2' <- get_sign_scores(calu.scov2.inf,
                                             gene_set = scov2_genes,
                                             slot = 'scale.data')

set.seed(333)
input <- FetchData(calu.scov2.inf, 
                   vars = unique(c(filter_degs$gene,
                                   'sarscov2')))
subsample <- sample(1:nrow(input), 3280)

## Random Forest
colnames(input) <- gsub('-', '_', colnames(input))
rforest <- randomForest(sarscov2 ~., data = input[subsample, ])
forest_pred <- predict(rforest, newdata = input[-subsample,])
pdf(paste0(path2_sc_deg, 'rforest_var_importance.pdf'))
plot(input$sarscov2[-subsample], forest_pred)
abline(1,1, col='red', lwd=3)
varImpPlot(rforest)
dev.off()

############################################################

signature.sel <- mutate(signature.sel, 
                        gene=rownames(signature.sel)) 
lapply(
        select(signature.sel, go_viral:infection_signature), 
        function(x) signature.sel$gene[x]
)

dir.create('data/22927382/')
## Downloading data
url_lorenz <- 'https://ndownloader.figshare.com/files/22927382'
lorenz_file <- 'data/Lorenz/covid_nbt_main.rds'
if ( ! file.exists(lorenz_file)) {
        download.file(url = url_lorenz, 
                      destfile = lorenz_file)
}

## Reading Lorenz data
lorenz <- readRDS('data/Lorenz/covid_nbt_main.rds')

## Permissivity 
## Reading permissivity signature
perm_signature <- read.table(
        'analysis/permissivity_signature_annotated.tsv',
        header = TRUE
)
head(perm_signature)
indxs <- 0.25*nrow(perm_signature) %>% as.integer()
indxs
signature <- perm_signature$gene[1:indxs]


## Scoring cells
## Cell scoring 
lorenz <- AddModuleScore(
        lorenz, 
        features = list(signature),
        name = 'permissivity_signature_score', 
        nbin = 200
)


pdf('figures/lorenz_scoring_by_celltype.pdf')
subset(lorenz, 
       severity %in% c('critical', 
                       'moderate'))@meta.data %>%
        ggplot(aes(x = celltype, 
                   y = permissivity_signature_score1,
                   fill = severity)) +
        scale_fill_manual(values = c('green', 'red')) +
        geom_boxplot() +
        coord_flip() +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        xlab('') + ggtitle('Lorenz')
dev.off()

pdf('figures/lorenz_scoring.pdf')
subset(lorenz, 
       severity %in% c('critical', 
                       'moderate'))@meta.data %>%
        ggplot(aes(x = severity, 
                   y = permissivity_signature_score1,
                   fill = severity)) +
        scale_fill_manual(values = c('green', 'red')) +
        geom_boxplot() +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        xlab('') + ggtitle('Lorenz')
dev.off()

pdf('figures/lorenz_scoring_boxplot_wt.pdf')
subset(lorenz, 
       severity %in% c('critical', 
                       'moderate'))@meta.data %>%
        mutate(severity = factor(severity,
                                 levels = c('moderate',
                                            'critical'))) %>%
        ggboxplot(x = "severity", 
                  y = "permissivity_signature_score1",
                  color = "severity", 
                  palette = "jco") +
        stat_compare_means(label.y = 0.8) +
        scale_color_manual(values = c('green3', 'red3')) +
        ggtitle('Lorenz')
dev.off()

#######################################################################
## Blinz data

blish <- readRDS('data/blish_covid.seu.rds')

## scoring cells
blish <- AddModuleScore(
        blish, 
        features = list(signature),
        name = 'permissivity_signature_score', 
        nbin = 50
)

pdf('figures/blish_scoring_by_celltype.pdf')
subset(blish, 
       Admission.level %in% c('Floor', 
                              'ICU'))@meta.data %>%
        ggplot(aes(x = cell.type.fine, 
                   y = permissivity_signature_score1,
                   fill = Admission.level)) +
        scale_fill_manual(values = c('green', 'red')) +
        geom_boxplot() +
        coord_flip() +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        xlab('') + ggtitle('Blish')
dev.off()

pdf('figures/blish_scoring.pdf')
subset(blish, 
       Admission.level %in% c('Floor', 
                              'ICU'))@meta.data %>%
        ggplot(aes(x = Admission.level, 
                   y = permissivity_signature_score1,
                   fill = Admission.level)) +
        scale_fill_manual(values = c('green', 'red')) +
        geom_boxplot() +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        xlab('') + ggtitle('Blish')
dev.off()

pdf('figures/blish_scoring_dotplot.pdf')
DotPlot(subset(blish, 
               Admission.level %in% c('Floor', 
                                      'ICU')),
        features = 'permissivity_signature_score1',
        group.by = 'Admission.level', 
        dot.scale = 20) +
        scale_color_viridis() + ggtitle('Blish')
dev.off()

pdf('figures/blish_scoring_boxplot_wt.pdf')
subset(blish, 
       Admission.level %in% c('Floor', 
                              'ICU'))@meta.data %>%
        mutate(Admission.level = factor(Admission.level,
                                        levels = c('Floor',
                                                   'ICU'))) %>%
        ggboxplot(x = "Admission.level", 
                  y = "permissivity_signature_score1",
                  color = "Admission.level", 
                  palette = "jco") +
        stat_compare_means(label.y = 0.11) +
        scale_color_manual(values = c('green3', 'red3')) +
        ggtitle('Blish')
dev.off()
