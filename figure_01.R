## Evaluation of the permissivity signature in mild vs severe
## cases
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

set.seed(333)

## Setting up paths
project_dir <- '~/sc/covid19-comorbidity/'
setwd(project_dir)

################################################################################
# Pie chart for all the top 25% most important genes in varios categories      #
################################################################################

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
