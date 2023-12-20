# CTC-smartseq3

Scripts required to replicate analysis of data for the following publication:


### File descriptions ###
environment.yml
- yaml file used to create conda environment in which the analysis was run

environment_exported.yml
- yaml file containing all versions of all packages in conda environment used for these analyses
- exported from environment creating using the environment.yml file (above)

installRdependencies.txt
- R packages installed on the R terminal

processCTC.R
- script used for analysis of scRNA-seq data

ensembleIDtoGeneNameMap.txt
- file mapping mitochondrial gene ID to gene name
- input for processCTC.R