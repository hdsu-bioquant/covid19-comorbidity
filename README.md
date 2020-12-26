# Analysis reproducibility

This repository provides all the scripts to reproduce the figures
of the manuscript: 'SPINT2 controls SARS-CoV-2 viral infection and 
is associated to disease severity'.

Custom R scripts are provided in the files *figure_01-02*, *03-04* and
*supp_Figures_01-03* to reproduce corresponding figures in the main
text and supplementary figures. 

The bash scripts *run_tobias* and *tobias_plot_tracks* shows how to
run the Transcription Factor footprinting analysis and plot
the genomic regions for SPINT2 and TMPRSS2 genomic loci in figure 1A,
respectively. Files to create the snakemake environment for running
tobias are provided in the envs folder and config files are given
in the configs directory.
