## Evaluation of the correlation of the colon organoids
source('https://raw.githubusercontent.com/KlugerLab/ALRA/master/alra.R')
source('https://raw.githubusercontent.com/KlugerLab/ALRA/master/alraSeurat2.R')
library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)
library(ggrepel)
library(scDoc)
library(scImpute)

#proyect_dir <- '/Users/carlosramirez/sc/serine_proteases_regulation/'
project_dir <- '/media/ag-cherrmann/cramirez/covid19-comorbidity/'
setwd(project_dir)

######################################################################
##                                                                  ##
##      Figure 4A. Correlation of SPRGs to viral load in patients   ##
##      nasopharynx samples from COVID-19 patients                  ##
##                                                                  ##
######################################################################

## Downloading data
url_lorenz <- 'https://ndownloader.figshare.com/files/22927382'
lorenz_file <- 'data/Lorenz/covid_nbt_main.rds'
if ( ! file.exists(lorenz_file)) {
    download.file(url = url_lorenz, 
                  destfile = lorenz_file)
}

## Reading Lorenz data
lorenz <- readRDS('data/Lorenz/covid_nbt_main.rds')

## Saving results to imput
lorenz@assays$RNA@counts %>%
    saveRDS(file = 'data/Lorenz/lorenz_raw_counts.rds')

lorenz.sp <- lapply(unique(lorenz@meta.data$sample), 
                    function(x) subset(lorenz, sample == x))
names(lorenz.sp) <- unique(lorenz@meta.data$sample)

## Data imputation
dir.create('data/tmp/')
for ( s in names(lorenz.sp)) {
    
    tmp_file <- paste0('data/tmp/', s, '.rds')
    saveRDS(lorenz.sp[s][[1]]@assays$RNA@counts,
            tmp_file)
    
    seurat <- readRDS(tmp_file)
    
    tmp_dir <- paste0('data/tmp/', s)
    dir.create(tmp_dir)
    
    scimpute(# full path to raw count matrix
        count_path = tmp_file, 
        infile = "rds",           # format of input file
        outfile = "rds",          # format of output file
        out_dir = tmp_dir,           # full path to output directory
        labeled = FALSE,          # cell type labels not available
        drop_thre = 0.5,          # threshold set on dropout probability
        Kcluster = length(unique(lorenz.sp[s][[1]]@meta.data$celltype)),             # 2 cell subpopulations
        ncores = 15) 
    
    file.remove(tmp_file)
    
}

## Integrating imputed results

## Getting output files names
imputed.names <- list.files('data/tmp', 
                            pattern = 'scimpute_count.rds', 
                            full.names = TRUE)

## Reading into R and transform to Seurat 
imputed.matrices <- sapply(imputed.names, readRDS)
lapply(imputed.matrices, function(x) x[1:3, 1:3])
imputed.seurat <- sapply(imputed.matrices, function(x)
    CreateSeuratObject(counts = x,
                       project = 'COVID19',
                       assay = 'RNA', 
                       min.cells = 1, 
                       min.features = 1)                 
)
imputed.seurat
rm(imputed.matrices)

## Integration
seurat.merged <- merge(x=imputed.seurat[1][[1]], 
                       y=imputed.seurat[2:length(imputed.seurat)])
seurat.merged@meta.data %>% head
seurat.merged.md <- seurat.merged@meta.data %>%
    add_rownames('id')

seurat.merged.md <- dplyr::select(seurat.merged.md, id)
## merge to lorenz metadata
lorenz.md <- lorenz@meta.data 
lorenz.md <- add_rownames(lorenz.md, var = 'id')
seurat.merged.md <- merge(lorenz.md, seurat.merged.md)
dim(seurat.merged.md)
dim(lorenz)
## Checking file names pairing
any(! colnames(seurat.merged) == seurat.merged.md$id )
## Correcting order
rownames(seurat.merged.md) <- seurat.merged.md$id
seurat.merged.md <- seurat.merged.md[colnames(seurat.merged), ]
any(! colnames(seurat.merged) == seurat.merged.md$id )
seurat.merged@meta.data <- seurat.merged.md

## Normalize data
seurat.merged <- NormalizeData(seurat.merged)
seurat.merged <- ScaleData(seurat.merged)

## Visualisation of the results

## Loading SPRGs
uniprot_SPRG.filtered <- readRDS('data/serine_proteases_filtered.rds')
SPRGs <- strsplit(uniprot_SPRG.filtered$Gene.names, split = ' |/') %>%
    unlist() %>% sort() %>% unique()
SPRGs

sprgs.df <- FetchData(subset(seurat.merged,
                             `SARS-CoV-2-p` > 0 | `SARS-CoV-2-m` > 0), 
                      vars = c(SPRGs, "SARS-CoV-2-p"))
sprgs.cor <- cor(sprgs.df)

pdf('figures/SPRGs_cor_viral_lorenz.pdf',
    height = 7, width = 7)
sprgs.cor %>%
    as.data.frame() %>%
    dplyr::select(`SARS-CoV-2-p`) %>%
    arrange(desc(abs(`SARS-CoV-2-p`))) %>%
    add_rownames(var = 'gene') %>%
    filter(gene != 'SARS-CoV-2-p') %>%
    mutate(rank=1:n()) %>%
    mutate(highlight=ifelse(rank < 20, 
                            TRUE,
                            FALSE)) %>%
    mutate(gene_label=ifelse(highlight ==TRUE,
                             gene, '')) %>%
    mutate(positive=ifelse(`SARS-CoV-2-p` >= 0,
                           'Posivitely correlated', 
                           'Negatively correlated')) %>%
    filter(!is.na(`SARS-CoV-2-p`)) %>%
    ggplot(aes(x=rank, 
               y=abs(`SARS-CoV-2-p`),
               label =gene_label,
               colour=positive)) +
    geom_point() +
    geom_text_repel(colour='black',
                    size=5) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    xlab('Rank') +
    ylab('Cor(SPINT2, viral P protein SARS-CoV2)') +
    scale_color_manual(values = c('red',
                                  'chartreuse4'))
dev.off()

######################################################################
##                                                                  ##
##      Figure 4B. ACE2, SPINT2 and TMPRSS2 gene expression in      ##
##      nasopharynx samples from COVID-19 patients                  ##
##                                                                  ##
######################################################################

## Mild vs Severe scRNASeq Lorenz et al, 2020 data
naso <- readRDS('~/Downloads/12436517/covid_nbt_main.rds')
naso

secretorybyseverity <- DotPlot(subset(naso, 
                                      celltype == 'Secretory' &
                                          severity != 'control'), 
                               group.by = 'severity', 
                               features = markers, 
                               dot.min = 0.1, 
                               col.min = 0, 
                               dot.scale = 35) +
    scale_colour_viridis() +
    theme_bw(base_size = 9) + 
    theme(panel.grid = element_blank(), 
          axis.text = element_text(colour="black", size=14), 
          axis.line = element_blank(), 
          text = element_text(family = 'Helvetica')) +
    xlab('') + ylab('') + coord_flip()

pdf('/Users/carlosramirez/sc/sars-cov2/figures/figure_02_00_mild_vs_severe_lorenz.pdf', 
    width = 15)
secretorybyseverity
dev.off()

######################################################################
##                                                                  ##
##      Figure 4C. ACE2, SPINT2 and TMPRSS2 genes expression        ##
##      in PBMC from COVID-19 patients                              ##
##                                                                  ##
######################################################################

## Mild vs severe scRNA-Seq Blish 
pbmc <- readRDS('~/Downloads/blish_covid.seu.rds')
pbmc <- subset(pbmc, Admission.level %in% c('ICU', 'Floor'))

pbmc.plot.subset <- DotPlot(subset(pbmc, cell.type %in% c('DC', 'pDC',
                                                          'SC & Eosinophil')), 
                            group.by = 'Admission.level', 
                            features = markers, 
                            dot.min = 0.1, 
                            col.min = 0, 
                            dot.scale = 50) +
    scale_colour_viridis() +
    theme_bw(base_size = 9) + 
    theme(panel.grid = element_blank(), 
          axis.text = element_text(colour="black", size=14), 
          axis.line = element_blank(), 
          text = element_text(family = 'Helvetica')) +
    xlab('') + ylab('') + 
    ggtitle('DCs + pDC + SC & Eosinophil') + coord_flip()
