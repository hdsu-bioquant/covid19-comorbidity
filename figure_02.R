## This script reproduce plots in figure 2 of the main 
## manuscript

## Dependencies
library(dplyr)
library(cmapR)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(ggrepel)
library(ggradar)


#######################################################################
##                                                                   ##                     
##      Figure 2A. Inference of the permissivity signature (i') =    ##
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

## Filtering induced genes upon infection
filter_degs <- filter(deg.markers, ( ! gene %in% infection.signature$gene) &
                              ( ! gene %in% h1299.infection.signature$gene) )

filter_degs %>% 
        saveRDS(paste0(path2_sc_deg, 'deg.calu_vs_h1299_mock.4hr.filtered.rds'))


##################################################################
##                                                              ##
##      Figure 1B. Ranking of the genes in the permissivity     ##
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

##################################################################################
##                                                                              ##
##      Figure 2B. Pie chart for all the top 25% most important genes           ##
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

#######################################################################
##                                                                   ##
##      Figure 2C. PEA of the top 25 % RF ranked genes in the        ##
##      permissivity signature                                       ##
##                                                                   ##
#######################################################################

enrRF = enrichr(genes = rf$Gene[rf$Weight > quantile(rf$Weight, 0.75)], 
                databases = c("GO_Biological_Process_2018")) #"KEGG_2019_Human", 
enrRF = reshape2::melt(lapply(enrRF, function(x) x[x$Adjusted.P.value < 0.05,][,c(1,2,8,4,9)]), 
                       id = c("Term", "Overlap", "Adjusted.P.value", "Combined.Score", "Genes"))                
enrRF = enrRF[,c(1,6,2,3,4,5)]
enrRF = enrRF[order(enrRF$Adjusted.P.value),]

tmp =  enrRF[1:10,c(1,4)]
tmp$Score = -log10(tmp$Adjusted.P.value)
tmp$Order = 2
tmp$Order[grep("viral", tmp$Term)] = 1
tmp$Order[grep("mRNA", tmp$Term)] = 3
tmp = tmp[order(tmp$Order, tmp$Adjusted.P.value),]
tmp$OrderVis = nrow(tmp):1

p2 = ggplot(dat = tmp, aes(x = reorder(Term, OrderVis), y = Score)) +  theme_bw(base_size = 9) +
        labs( x = "", y = "-log10 Adjusted Pvalue") +
        geom_bar(stat="identity", fill = c(rep("#9e9ac8", 2), rep("#756bb1", 5), rep("#cbc9e2", 3))) +
        geom_hline(yintercept = -log10(0.01)) +
        coord_flip() +
        theme(panel.grid = element_blank(),
              axis.text = element_text(colour="black"), 
              axis.line = element_blank())
rm(tmp)

pdf('figures/enrichr.pdf')
p2
dev.off()

######################################################################
##                                                                  ##
##      Figure 2D. Scoring permissivity signature in HCL cells      ##
##                                                                  ##
######################################################################

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

######################################################################
##                                                                  ##
##      Figure 2E. Correlations of top % 25 RF ranked genes         ##
##      to SARS-CoV2 viral load translation rates in Caco-2         ##
##      cell lines                                                  ##
######################################################################

protdat = rio::import("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2332-7/MediaObjects/41586_2020_2332_MOESM2_ESM.xlsx")
table(protdat$`Species Names01`)

# Viral proteins 
vp = protdat[protdat$`Species Names01` == "Wuhan seafood market pneumonia virus OX=2697049",]
rownames(vp) = vp$Accession
vp = t(vp[,16:27])

cor(vp, method="spearman")
pheatmap(cor(vp, method="spearman"))
#PODTC7 is an outlier
#PODTC2 and PODTC8 cluster together and
#PODTS2 and PODTC9 cluster together
#Selecting one viral protein (Spike and Nucleoprotein) from each cluster

# Permissive genes from RF analysis in the proteome data
pg = protdat[protdat$`Gene Symbol01` %in% rf$Gene[rf$Weight > quantile(rf$Weight, 0.75)],]
rownames(pg) = pg$`Gene Symbol01`
pg = pg[,16:27]
pg = pg[rowSums(pg) > 0,]

length(rf$Gene[rf$Weight > quantile(rf$Weight, 0.75)])
#[1] 118
nrow(pg)
#[1] 56

# Correlation between Viral proteins and Permissive genes
corRes = data.frame(matrix(ncol = ncol(vp), nrow= nrow(pg)))
colnames(corRes) = colnames(vp)

for(i in 1:nrow(corRes))
{
        for(j in 1: ncol(vp))
        {
                corRes[i,j] = cor(as.numeric(vp[,j]), 
                                  as.numeric(pg[i,]), 
                                  method="spearman", 
                                  use="pairwise.complete.obs")
        }
}

corRes$Gene = rownames(pg)
corRes$Label = "No"
corRes$Label[abs(corRes$P0DTC2) > 0.5 & abs(corRes$P0DTC9) > 0.5] = "Yes"
corRes = merge(x = corRes, y = gs, by = "Gene", all.x = T)

corRes$group[is.na(corRes$group)] = "Other"
corRes$color[corRes$group == "Other"] = "black"

p4 = ggplot(corRes, aes(x = P0DTC2, y = P0DTC9)) + 
        labs(x = "P0DTC2 | Spike glycoprotein SARS-CoV2", y = "P0DTC9 | Nucleoprotein SARS-CoV2") +
        geom_point(color = "#756bb1") + theme_bw(base_size = 9) +
        geom_vline(xintercept = c(-0.5,0.5)) + geom_hline(yintercept = c(-0.5,0.5)) +
        theme(panel.grid = element_blank(), axis.text = element_text(colour="black"), axis.line = element_blank()) + 
        geom_text_repel(data = subset(corRes, Label == "Yes"), aes(label = Gene), 
                        colour = subset(corRes, Label == "Yes")$color, size = 2)

