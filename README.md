# FAST Genomics-compatible single-cell RNA-seq preprocessing and clustering
This container executes an R script that uses functionality provided by the Seurat package (http://satijalab.org/seurat/) to perform initial quality control, selection of variable genes, dimensionality reduction via PCA and tSNA clustering of single-cell RNA-seq data.
## Data In-/Output
All data read and written by the analysis script comes from / goes to the folder that has to be mounted as '/data'. 
Single-cell RNA-seq data has to be provided in the sparse matrix FAST Genomics format in a file called 'expressions_entrez.tsv'.
Clustering output is provided in the comma-separated file 'clustering_results.csv'.
