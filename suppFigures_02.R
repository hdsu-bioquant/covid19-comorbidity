## Scripts to reproduce supplementary figure 2 of the manuscript
## SPINT2 controls SARS-CoV-2 viral infection and is associated to disease severity

library(Seurat)
library(viridis)
library(ggplot2)
library(tidyverse)

markers <- c('SPINT2',
             'ACE2', 'TMPRSS2')

path2project <- '/Users/carlosramirez/sc/sars-cov2'
setwd(path2project)

#######################################################################
##                                                                   ##
##      Supp figure 2A. ACE2, SPINT2 and TMPRSS2 gene expression     ##
##      in nasopharynx samples from COVID-19 patients                ##
##                                                                   ##
#######################################################################

naso <- readRDS('data/22927382/covid_nbt_main.rds')

DotPlot(naso, 
        group.by = 'celltype', 
        features = c('ACE2', 'TMPRSS2', 'SPINT2')) +
        scale_colour_viridis() +
        theme_bw(base_size = 9) + 
        theme(panel.grid = element_blank(), 
              axis.text = element_text(colour="black",
                                       size=14), 
              axis.line = element_blank()) +
        xlab('') + ylab('')

#######################################################################
##                                                                   ##
##      Supp Figure 2B. ACE2, SPINT2 and TMPRSS2 gene expression     ##
##      in PBMC cells from COVID-19 patients                         ##
##                                                                   ##
#######################################################################

## Wilk AJ. et al, 2020. Expression of SPINT2
pbmc <- readRDS('data/blish_covid.seu.rds')
pbmc <- subset(pbmc, Admission.level %in% c('ICU', 'Floor'))

## ACE2, TMPRSS2, and SPINT2 expression in PBMC cells
pdf('~/sc/sars-cov2/figures/spint2_blish.pdf')
DotPlot(pbmc, 
        group.by = 'cell.type', 
        features = sarscov2_markers, 
        dot.min = 0.1, 
        col.min = 0, 
        dot.scale = 30) +
        scale_colour_viridis() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        xlab('') + ylab('')
dev.off()

#######################################################################
##                                                                   ##
##      Supp figure 2C. Correlation of ACE2, SPINT2 and TMPRSS2      ##
##      to viral gene expression in COVID-19 patients                ##
##                                                                   ##
#######################################################################

protdat = rio::import("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2332-7/MediaObjects/41586_2020_2332_MOESM2_ESM.xlsx")
table(protdat$`Species Names01`)

# Viral proteins 
vp = protdat[protdat$`Species Names01` == "Wuhan seafood market pneumonia virus OX=2697049",]
rownames(vp) = vp$Accession
vp = t(vp[,16:27])

cor(vp, method="spearman")
pheatmap(cor(vp, method="spearman"))

# Permissive genes from RF analysis in the proteome data
pg = protdat[protdat$`Gene Symbol01` %in% rf$Gene[rf$Weight > quantile(rf$Weight, 0.75)],]
rownames(pg) = pg$`Gene Symbol01`
pg = pg[,16:27]
pg = pg[rowSums(pg) > 0,]

## Reading results from RF
rf <- read.table('data/variable_importance_ranked_genes.tsv',
                 header = TRUE)

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


p_up = ggarrange(p1, p2, ncol=2, widths = c(0.45,1))
p_dw = ggarrange(p3, p4, ncol=2, widths = c(0.8,1))
ggarrange(p_up, p_dw, ncol = 1) %>% 
        ggexport(filename = "figures/Analysis_of_genes_associated_with_viral_load.pdf", width = 8, height = 5)

corRes = corRes[!duplicated(corRes$Gene),]
tmp = corRes[,2:6]
rownames(tmp) = corRes$Gene
tmpanno = data.frame(corRes[,c("Label", "group", "color")], stringsAsFactors = F)
rownames(tmpanno) =  corRes$Gene
tmpcol = list(Label = setNames(c("black", "grey80"), c("Yes", "No")),
              group = sapply(split(tmpanno$color, tmpanno$group), unique))
tmpcol$group[tmpcol$group == "black"]="grey80" 
tmpanno = tmpanno[,1:2]

pheatmap(tmp, clustering_method = "ward.D2",
         breaks = seq(-1,1,0.25),
         color = rev(brewer.pal(10, "RdYlBu")),
         annotation_row = tmpanno, treeheight_row = 0, treeheight_col = 0,
         cutree_rows = 3, annotation_colors  = tmpcol,
         border_color = NA, fontsize = 6, 
         filename = paste0(path, "results/Analysis_of_genes_associated_with_viral_load_all_genes.pdf"),
         width = 3.5, height = 5)


