## GTEx exploration of the gene expression profiles of serine 
## proteases
library(cmapR)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(ggrepel)
library(DESeq2)

## Setting paths
setwd('/Users/carlosramirez/sc/serine_proteases_regulation/')

url_file <- 'https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz'
file_path <- 'data/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct'

## Downloading 
if ( ! file.exists(file_path) ){
        download.file(url = url_file,
                      destfile = file_path)
}

#gtex <- parse_gctx(file_path)
#dim(gtex@mat)
#gtex@mat[1:5, 1:5]
#saveRDS(gtex@mat, 'data/gtex/GTEx_tpm.mtx')
#rm(gtex)

#gtex.mtx <- readRDS('data/gtex/GTEx_tpm.mtx')
#gtex.mtx[1:5, 1:5]

#############################################################
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

head(gene_anns)

## Cleanning gene names
#gtex.mtx <- gtex.mtx[rownames(gtex.mtx) %in% 
#                             rownames(gene_anns), ] 
#gtex.mtx <- gtex.mtx[rownames(gene_anns), ]
#rownames(gtex.mtx) <- gene_anns$gene_name
#gtex.mtx[1:5, 1:5]

#saveRDS(gtex.mtx, 'data/gtex/GTEx_tpm.mtx')
###########################################################
gtex.mtx <- readRDS('data/gtex/GTEx_tpm.mtx')

## Sample annotations
gtex.ann <- read.table(
        'https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt',
        header = TRUE, sep = '\t',
        quote = ''
)
dim(gtex.ann)
gtex.ann[1:5, 1:5]

#############################################################
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
#rm(gtex.mtx)

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


########################################################
## Enterocytes analysis
ent.gtex <- dplyr::filter(ge.ann, SMTS == "Small Intestine")

sprs.cor <- cor(ent.gtex[, colnames(ent.gtex) %in% SPRGs], 
                method = 'spearman')
diag(sprs.cor) <- 0

#######################################################
##                                                   ##
##      Heatmap of highly correlated SPRGs           ##
##                                                   ##
#######################################################
## SPRGs
pdf('figures/SPRG_correlation_ileum_GTEx.pdf',
    width = 16, height = 24)
pheatmap::pheatmap(sprs.cor, color = viridis(20))
dev.off()

###########################################################
## SPINT2 correlation

## TFs 
tf.cor <- cor(ent.gtex[, colnames(ent.gtex) %in% 
                               c(tfs, 
                                 'SPINT2', 
                                 'TMPRSS2')], 
              method = 'spearman')
diag(tf.cor) <- 0

spint2.cor.df <- data.frame(spint2_correlation=tf.cor[, 'SPINT2'])
head(spint2.cor.df)

mean.df <- data.frame(mean_val=sapply(ent.gtex, 
                                      mean, 
                                      na.rm=TRUE))
head(mean.df)

vars <- intersect(rownames(spint2.cor.df), rownames(mean.df))
spint2.cor.df <- as.data.frame(spint2.cor.df[vars,])
rownames(spint2.cor.df) <- vars
mean.df <- mean.df[vars,]
## check rownames are equal
any(rownames(spint2.cor.df)!=rownames(mean.df))
df <- cbind(spint2.cor.df, mean.df)
rownames(df) <- rownames(spint2.cor.df)
colnames(df) <- c('correlation', 'mean')
head(df)

## Thresholds
cor_thr <- 0.8
mean_thr <- 2.5
xlims <- 1.05

df <- df %>%
        as.data.frame() %>%
        mutate(gene=rownames(df)) %>%
        mutate(highlight=ifelse(log(mean+1) > mean_thr & 
                                        ( correlation > 0.8  |
                                                  correlation < - 0.85),
                                TRUE, FALSE)) %>%
        mutate(gene_label=ifelse(highlight==TRUE, 
                                 gene, ''))

cor_genes.ileum <- filter(df, highlight == TRUE &
                                  correlation > 0 ) %>%
                        select(gene) %>%
                        unlist() %>%
                        unname()
                
        
pdf('figures/spint2_tfs_correlation_ileum_gtex.pdf',
    width = 9, height = 7)
spint2.plot <- df %>%
        ggplot(aes(x=correlation, y=log(mean+1),
                   colour=highlight, 
                   label=gene_label)) +
        geom_point() +
        geom_vline(xintercept = c(-cor_thr, cor_thr),
                   linetype = 'dashed',
                   colour='gray5') +
        geom_hline(yintercept = mean_thr,
                   linetype = 'dashed',
                   colour='gray5') +
        geom_text_repel() +
        scale_color_manual(values = c('black', 'red')) +
        theme_bw() +
        theme(legend.position = 'none',
              panel.grid = element_blank()) +
        xlim(-xlims, xlims) +
        xlab('Spearman Rho values') +
        ylab('Log(mean TPM + 1)')
spint2.plot
dev.off()

##############################################################
## TMPRSS2
tmprss2.cor.df <- data.frame(tmprss2.df_correlation=tf.cor[, 'TMPRSS2'])
head(tmprss2.cor.df)

mean.df <- data.frame(mean_val=sapply(ent.gtex, 
                                      mean, 
                                      na.rm=TRUE))
head(mean.df)

vars <- intersect(rownames(tmprss2.cor.df), rownames(mean.df))
tmprss2.cor.df <- as.data.frame(tmprss2.cor.df[vars,])
rownames(tmprss2.cor.df) <- vars
mean.df <- mean.df[vars,]
## check rownames are equal
any(rownames(tmprss2.cor.df)!=rownames(mean.df))
tmprss2.df <- cbind(tmprss2.cor.df, mean.df)
rownames(tmprss2.df) <- rownames(tmprss2.cor.df)
colnames(tmprss2.df) <- c('correlation', 'mean')
head(tmprss2.df)

tmprss2.df <- tmprss2.df %>%
        as.data.frame() %>%
        mutate(gene=rownames(tmprss2.df)) %>%
        mutate(highlight=ifelse(log(mean+1) > mean_thr & 
                                        ( correlation > 0.8  |
                                                  correlation < - 0.85),
                                TRUE, FALSE)) %>%
        mutate(gene_label=ifelse(highlight==TRUE, 
                                 gene, '')) %>%
        filter( ! gene %in% c('SPINT2', 'TMPRSS2'))
head(tmprss2.df)

cor_tmprss2_genes.ileum <- filter(tmprss2.df, 
                                  highlight == TRUE &
                                  correlation > 0 ) %>%
        select(gene) %>%
        unlist() %>%
        unname()


cor_thr <- 0.8
mean_thr <- 2.5
xlims <- 1.05
pdf('figures/tmprss2_tfs_correlation_ileum_gtex.pdf',
    width = 9, height = 7)
tmprss2.plot <- tmprss2.df %>%
        ggplot(aes(x=correlation, y=log(mean+1),
                   colour=highlight, 
                   label=gene_label)) +
        geom_point() +
        geom_vline(xintercept = c(-cor_thr, cor_thr),
                   linetype = 'dashed',
                   colour='gray5') +
        geom_hline(yintercept = mean_thr,
                   linetype = 'dashed',
                   colour='gray5') +
        geom_text_repel() +
        scale_color_manual(values = c('black', 'red')) +
        theme_bw() +
        theme(legend.position = 'none',
              panel.grid = element_blank()) +
        xlim(-xlims, xlims) +
        xlab('Spearman Rho values') +
        ylab('Log(mean TPM + 1)')
tmprss2.plot
dev.off()

###############################################################

plot <- df %>%
        mutate(highlight = ( highlight == TRUE & 
                                     gene_label %in% cor_tmprss2_genes.ileum )) %>%
        ggplot(aes(x=correlation, y=log(mean+1),
                   colour=highlight, 
                   label=gene_label)) +
        geom_point() +
        geom_vline(xintercept = c(-cor_thr, cor_thr),
                   linetype = 'dashed',
                   colour='gray5') +
        geom_hline(yintercept = mean_thr,
                   linetype = 'dashed',
                   colour='gray5') +
        geom_text_repel() +
        scale_color_manual(values = c('black', 'red')) +
        theme_bw() +
        theme(legend.position = 'none',
              panel.grid = element_blank()) +
        xlim(-xlims, xlims) +
        xlab('Spearman Rho values') +
        ylab('Log(mean TPM + 1)')
plot

pdf('figures/tfs_correlation_ileum_gtex.pdf',
    width = 15)
gridExtra::grid.arrange(
        spint2.plot, 
        tmprss2.plot, 
        ncol = 2
)
plot
dev.off()

