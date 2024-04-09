# Gilia-QTL-data
Scripts and data (phenotypes and genotypes) for running QTL analysis in Gilia

All scripts and data files should be in the same folder/working directory when running

## Scripts:

QTL_dataprep.R—Reads in the "cross" object required for R/qtl analyses, formats trait names appropriately, and creates the interval table for significant QTL

CIM_interval_plot.R—Shows marker density along each chromosome, positions of each QTL, and trait names above QTL for reference.

2x2_QTL_graphs.R—Generates the QTL traces present in the manuscript in 2x2 format.

Histogram and Correlation Figures.R—Generates the histogram figure and the correlation figure from the manuscript.

Trait_means_table.R—Generates the phenotypic data tables present in the manuscript.

Normality tests.R—Tests for normality using the Shapiro-Wilks test (several other methods were also tested in the script).
