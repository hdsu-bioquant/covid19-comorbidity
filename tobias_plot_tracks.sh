## Command to generate Figure 1b. TFs bound to
## SPINT2 and TMPRSS2 genomic loci 
TOBIAS PlotTracks \
       --bigwigs analysis/TOBIAS/hio_footprints.bw \
       --regions analysis/gene_coordinates.tsv \
       --sites analysis/binding_sites.tsv \
       --gtf data/gencode.v35.annotation.gtf \
       --colors red darkblue
