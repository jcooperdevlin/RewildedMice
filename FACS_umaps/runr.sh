#!/bin/bash

#SBATCH -t 10:00:00 
#SBATCH -c 1	
#SBATCH --mem-per-cpu 50GB
#SBATCH --job-name runr
#SBATCH -o stdout.txt
#SBATCH -e error.txt

module load r/3.4.2

#Rscript flow_combiner2.20.R
#Rscript plot_tsne_cluster.R
#Rscript plot_all.R
#Rscript plot_combo.R

Rscript Step3_plot_FACS.R
