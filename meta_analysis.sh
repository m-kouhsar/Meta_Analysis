#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=10:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --output=%x.%j.out
#########################################################################
# Input arguments:
#   
#   Input_Dir: A directory contains Input files 
#   
#   Input_csv_file: A csv (comma delimited) files must contains the following columns (column names are NOT case sensitive).
#             Study_name: Name of study
#             DEG_file: Full path of DEG results file for each study in the Input_Dir.
#             Gene_col: Column name for gene IDs in DEG results (default is "Gene")
#             SE_col: column name for standard error values in DEG results ("Default is "SE")
#             ES_col: column name for effect size values in DEG results (Default is "ES")
#   
#   min_gene_appearance: minimum number times each gene appears across the studies (Default is number of studies or number of rows in Input_csv_file)

# For each study, a DEG results file is required. It must be a CSV file (comma delimited) contains at least 3 columns for genes, standard errors and estimate values.

###################################################################################

Input_Dir="./DEG_results/"
Input_csv_file="meta_input_exceRpt.Br0.Br2.csv"
min_gene_appearance=5
OutPrefix="./Results/meta_analysis/SmallRNAProject"

ScriptDir="/Morteza/github/Meta_Analysis/"
###################################################################################

Rscript ${ScriptDir}/meta_analysis.R "$Input_Dir"  "$Input_csv_file"  "$min_gene_appearance"  "$OutPrefix" 