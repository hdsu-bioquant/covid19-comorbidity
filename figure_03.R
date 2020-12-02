## Evaluation of the correlation of the colon organoids
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(viridis)
library(dorothea)
library(viper)

proyect_dir <- '/Users/carlosramirez/sc/serine_proteases_regulation/'
setwd(proyect_dir)

###############################################################
##                                                           ##
##      Figure 4c. TF activities in scRNA-Seq using          ##
##      Ileum organoid data.                                 ##
##                                                           ##
###############################################################


## Getting SPRG
uniprot_SPRG.filtered <- readRDS(
    'data/serine_proteases_filtered.rds'
)
SPRGs <- strsplit(uniprot_SPRG.filtered$Gene.names, 
                  split = ' |/') %>%
        unlist() %>% sort() %>% unique()

load('~/sc/organoids/data/COVID19_July.rda')

## Loading TF activities
TF_act <- read.csv(
        '~/sc/organoids/analysis/scenic/ileum_corrected/ileum_aucell.csv'
)
head(TF_act)
dim(TF_act)

## Getting serine protease gene expression
DefaultAssay(Illeum_H_T) <- 'RNA'
SPRG.df <- FetchData(Illeum_H_T, vars = SPRGs)
dim(SPRG.df)

## Checking rownames pairing
any(rownames(SPRG.df) != TF_act$Cell)

## Merging results
SPRG.df <- SPRG.df[, sapply(SPRG.df, sum) > 0]
ileum.cor <- cor(SPRG.df, select(TF_act, -Cell), 
                 method = 'spearman')
ileum.cor[1:5, 1:5]
colnames(ileum.cor) <- gsub('\\.', '', 
                            colnames(ileum.cor))
dim(ileum.cor)

## Calculating mean TF activities
TFA_mean <- apply(select(TF_act, -Cell), 2, mean)
names(TFA_mean) <-  gsub('\\.', '', names(TFA_mean))

## tfa correlation vs gene expression in ileum
spint2_cor_tfa <- ileum.cor['SPINT2', ]

vars <- intersect(names(TFA_mean),
                  names(spint2_cor_tfa))


ileum.tfa.df <- data.frame(
    spint2_cor_tfa = spint2_cor_tfa[vars],
    mean_TF_Activity = tfs.mean_exp[vars]
)
dim(ileum.tfa.df)

## Adding TMPRSS2 coregulation
tmprss2_cor_tfa <- ileum.cor['TMPRSS2', ]
ileum.tfa.df$'tmprss2_cor_tfa' <- tmprss2_cor_tfa[vars]

st14_cor_tfa <- ileum.cor['ST14', ]
ileum.tfa.df$'st14_cor_tfa' <- st14_cor_tfa[vars]

spint1_cor_tfa <- ileum.cor['SPINT1', ]
ileum.tfa.df$'spint1_cor_tfa' <- spint1_cor_tfa[vars]



cor_th <- 0
mean_th <- 0.5
pdf('figures/tfa_cor_spint2_tmprss2_ileum_org.pdf',
    height = 7, width = 8)
ileum.tfa.df %>%
        mutate(gene = rownames(ileum.tfa.df)) %>%
        mutate(highlight = ifelse( mean_TF_Activity > mean_th,
                                  TRUE, FALSE)) %>%
        mutate(gene_label = ifelse(highlight == TRUE,
                                   gene, '')) %>%
        ggplot(aes(x=spint2_cor_tfa, 
                   y=tmprss2_cor_tfa,
                   colour=mean_TF_Activity,
                   label = gene_label)) +
                geom_point(aes(size=mean_TF_Activity)) +
                geom_text_repel(colour='black') +
                theme_bw() +
                scale_color_viridis() +
                theme(panel.grid = element_blank()) +
                geom_hline(yintercept = cor_th, 
                           linetype='dashed') +
                geom_vline(xintercept = cor_th,
                           linetype = 'dashed') +
                ylab('Cor(TF activity, TMPRSS2)') +
                xlab('Cor(TF activity, SPINT2)') +
                ggtitle('Ileum organoid')
dev.off()


pdf('figures/tfa_cor_spint1_st14_ileum_org.pdf',
    height = 7, width = 8)
ileum.tfa.df %>%
    mutate(gene = rownames(ileum.tfa.df)) %>%
    mutate(highlight = ifelse( mean_TF_Activity > mean_th,
                               TRUE, FALSE)) %>%
    mutate(gene_label = ifelse(highlight == TRUE,
                               gene, '')) %>%
    ggplot(aes(x=spint1_cor_tfa, 
               y=st14_cor_tfa,
               colour=mean_TF_Activity,
               label = gene_label)) +
    geom_point(aes(size=mean_TF_Activity)) +
    geom_text_repel(colour='black') +
    theme_bw() +
    scale_color_viridis() +
    theme(panel.grid = element_blank()) +
    geom_hline(yintercept = cor_th, 
               linetype='dashed') +
    geom_vline(xintercept = cor_th,
               linetype = 'dashed')+
    ylab('Cor(TF activity, ST14)') +
    xlab('Cor(TF activity, SPINT1)') +
    ggtitle('Ileum organoid')
dev.off()


###############################################################
##                                                           ##
##      Figure 4d. TF Binding to SPINT2 and TMRPSS2          ##
##      as observed in ATAC-Seq data from Human Intes-       ##
##      tinal organoids                                      ##
##                                                           ##
###############################################################

## Loading TOBIAS results
files <- list.files('~/sc/covid19-comorbidity/analysis/TOBIAS/tfbs_HIO/', 
                    full.names = TRUE)
tfbs.list <- lapply(files,
                    function(x)
                        read.table(x, header = TRUE))

## Merging TFBS
tfbs.df  <- do.call(rbind, tfbs.list)

## Getting common TF binding to SPINT2 and TMPRSS2
spint2.tfb <- filter(tfbs.df, gene_name == 'SPINT2') %>%
    dplyr::select(TFBS_name) %>%
    unlist() %>% unique()
tmprss2.tfb <- filter(tfbs.df, gene_name == 'TMPRSS2') %>%
    dplyr::select(TFBS_name) %>%
    unlist() %>% unique()
tfbs.common <- intersect(spint2.tfb, tmprss2.tfb)
tfbs.common.df <- filter(tfbs.df,
                         gene_name %in% c('SPINT2', 
                                          'TMPRSS2')) %>%
    filter(TFBS_name %in% tfbs.common) 

## Filtering using most correlated genes in scRNA-Seq data
## from Figure 4c results
pos.cor <- ileum.tfa.df %>%
            filter(spint2_cor_tfa > 0 & 
                       tmprss2_cor_tfa > 0) %>%
                rownames()
pos.cor

tfbs.common.df %>%
    mutate(TFBS_name = gsub('_.*', '', TFBS_name)) %>%
    filter(TFBS_name %in% pos.cor) %>%
    group_by(TFBS_name, gene_name) %>%
    summarise(max=max(TFBS_score)) %>%
    mutate(TFBS_score=max) %>%
    mutate(TFBS_score_2=TFBS_score) %>%
    dplyr::select(TFBS_name, gene_name, TFBS_score, TFBS_score_2) %>%
    write.table('~/sc/covid19-comorbidity/analysis/network_tobias.tsv',
                sep = '\t',
                quote = FALSE,
                row.names = FALSE)

tfs <- gsub('_.*', '', tfbs.common.df$TFBS_name)
intersection <- intersect(pos.cor, 
                          unique(tfs))

cor_th <- 0
mean_th <- 0.5
pdf('figures/tfa_cor_spint2_tmprss2_ileum_org_v2.pdf',
    height = 7, width = 8)
ileum.tfa.df %>%
    mutate(gene = rownames(ileum.tfa.df)) %>%
    mutate(highlight = ifelse( gene %in% intersection,
                               TRUE, FALSE)) %>%
    mutate(gene_label = ifelse(highlight == TRUE,
                               gene, '')) %>%
    ggplot(aes(x=spint2_cor_tfa, 
               y=tmprss2_cor_tfa,
               colour=mean_TF_Activity,
               label = gene_label)) +
    geom_point(aes(size=mean_TF_Activity)) +
    geom_text_repel(colour='black') +
    theme_bw() +
    scale_color_viridis() +
    theme(panel.grid = element_blank()) +
    geom_hline(yintercept = cor_th, 
               linetype='dashed') +
    geom_vline(xintercept = cor_th,
               linetype = 'dashed') +
    ylab('Cor(TF activity, TMPRSS2)') +
    xlab('Cor(TF activity, SPINT2)') +
    ggtitle('Ileum organoid')
dev.off()

###########################################################
## Same analysis for ST14 and SPINT1

## Getting common TF binding to SPINT2 and TMPRSS2
spint1.tfb <- filter(tfbs.df, gene_name == 'SPINT1') %>%
    dplyr::select(TFBS_name) %>%
    unlist() %>% unique()
st14.tfb <- filter(tfbs.df, gene_name == 'ST14') %>%
    dplyr::select(TFBS_name) %>%
    unlist() %>% unique()
tfbs.common <- intersect(spint1.tfb, st14.tfb)
tfbs.common.df <- filter(tfbs.df,
                         gene_name %in% c('SPINT1', 
                                          'ST14')) %>%
    filter(TFBS_name %in% tfbs.common) 

## Filtering using most correlated genes in scRNA-Seq data
## from Figure 4c results
pos.cor <- ileum.tfa.df %>%
    filter(spint1_cor_tfa > 0 | 
               st14_cor_tfa > 0) %>%
    rownames()
pos.cor


tfbs.common.df %>%
    mutate(TFBS_name = gsub('_.*', '', TFBS_name)) %>%
    filter(TFBS_name %in% pos.cor) %>%
    group_by(TFBS_name, gene_name) %>%
    summarise(max=max(TFBS_score)) %>%
    mutate(TFBS_score=max) %>%
    mutate(TFBS_score_2=TFBS_score) %>%
    dplyr::select(TFBS_name, gene_name, TFBS_score, TFBS_score_2) %>%
    write.table('~/sc/covid19-comorbidity/analysis/network_tobias_st14.tsv',
                sep = '\t',
                quote = FALSE,
                row.names = FALSE)

tfs <- gsub('_.*', '', tfbs.common.df$TFBS_name)
intersection <- intersect(pos.cor, 
                          unique(tfs))
intersection

ileum.tfa.df %>%
    mutate(gene = rownames(ileum.tfa.df)) %>%
    mutate(highlight = ifelse( gene %in% intersection,
                               TRUE, FALSE)) %>%
    mutate(gene_label = ifelse(highlight == TRUE,
                               gene, '')) %>%
    ggplot(aes(x=spint1_cor_tfa, 
               y=st14_cor_tfa,
               colour=mean_TF_Activity,
               label = gene_label)) +
    geom_point(aes(size=mean_TF_Activity)) +
    geom_text_repel(colour='black') +
    theme_bw() +
    scale_color_viridis() +
    theme(panel.grid = element_blank()) +
    geom_hline(yintercept = cor_th, 
               linetype='dashed') +
    geom_vline(xintercept = cor_th,
               linetype = 'dashed') +
    ylab('Cor(TF activity, TMPRSS2)') +
    xlab('Cor(TF activity, SPINT2)') +
    ggtitle('Ileum organoid')
