## Evaluation of the permissivity signature in mild vs severe
## cases
library(Seurat)
library(dplyr)
library(ggplot2)

set.seed(333)

## Setting up paths
project_dir <- '~/sc/covid19-comorbidity/'
setwd(project_dir)

## Downloading data
url_lorenz <- 'https://ndownloader.figshare.com/files/22927382'
if ( ! dir.exists('data/22927382')) {
        download.file(url = url_lorenz, 
                      destfile = 'data/22927382/')
}

## Reading Lorenz data
lorenz <- readRDS('data/22927382/covid_nbt_main.rds')

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


## subsampling
lorenz.ss <- subset(lorenz, 
                    cells = sample(Cells(lorenz), 
                                   3000))
lorenz.ss$severity %>% table()
lorenz.ss$celltype %>% table()
rm(lorenz)

## Scoring cells
## Cell scoring 
lorenz.ss <- AddModuleScore(
        lorenz.ss, 
        features = list(signature),
        name = 'permissivity_signature_score', 
        nbin = 50
)

## Saving results
saveRDS(lorenz.ss@meta.data, 'analysis/lorenz_scoring.rds')
write.table(lorenz.ss@meta.data, 
            'analysis/lorenz.ss_scoring.tsv', 
            sep='\t')

subset(lorenz.ss, 
       severity %in% c('critical', 
                       'moderate'))@meta.data %>%
        ggplot(aes(x = celltype, 
                   y = permissivity_signature_score1,
                   fill = severity)) +
                geom_boxplot() +
                coord_flip() +
                theme_bw() +
                theme(panel.grid = element_blank()) +
                xlab('')

