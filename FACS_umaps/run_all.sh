#!/bin/bash

#SBATCH -t 12:00:00 
#SBATCH -c 1
#SBATCH --mem-per-cpu 50GB
#SBATCH --job-name fast_tsne
#SBATCH -o tsne_umap_stdout.txt
#SBATCH -e tsne_umap_error.txt

module load anaconda3/cpu/5.2.0
source activate /gpfs/data/ruggleslab/home/devlij03/my_env
module unload anaconda3/cpu/5.2.0
module load fftw3/openmpi/gcc/64/3.3.7 
module load r/3.4.2

## Organize Inputs
#Rscript Step1_combine_FACS.R

#python fast_tsne.py inputs/MLN_df.txt inputs/ftsne_combo_MLN.txt inputs/umap_combo_MLN.txt
#python fast_tsne.py inputs/MLN_myeloid_df.txt inputs/ftsne_combo_MLN_myeloid.txt inputs/umap_combo_MLN_myeloid.txt
#python fast_tsne.py inputs/Blood_df.txt inputs/ftsne_combo_Blood.txt inputs/umap_combo_Blood.txt

#Rscript Step2_subset_FACS.R


#python fast_umap.py inputs/CD4_Blood_df.txt inputs/CD4_umap_combo_Blood.txt
#python fast_umap.py inputs/CD8_Blood_df.txt inputs/CD8_umap_combo_Blood.txt
#python fast_umap.py inputs/CD19_Blood_df.txt inputs/CD19_umap_combo_Blood.txt

python fast_umap.py inputs/CD19_MLN_df.txt inputs/CD19_umap_combo_MLN.txt
#python fast_umap.py inputs/CD4_MLN_df.txt inputs/CD4_umap_combo_MLN.txt
#python fast_umap.py inputs/CD8_MLN_df.txt inputs/CD8_umap_combo_MLN.txt
