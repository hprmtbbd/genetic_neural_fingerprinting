# Overview
Custom analysis code for the manuscript entitled "Genetic fingerprinting with heritable phenotypes of the resting-state brain network topology," published in _Communications Biology_. All scripts were written in MATLAB R2020a or R2023b.

# Scripts
The mfile_h folder contains scripts for performing statistical analyses on the heritability values, generating figures and tables of the results, and implementing the genetic fingerprinting machine learning algorithm.
- fMRI_MEG_GGM_LGM_h2_analysis.m (written in MATLAB 2020a) will generate figures and tables of the heritability values for the global and local graph measures in the MEG and fMRI datasets. This includes Fig. 2-5 in the main text and Table S2-S3, Table S8-S9, and Fig. S1-S8 in the Supplementary Information.
- plot_heritability_pipeline_figures.m (written in MATLAB 2023b) will generate figures for comparing the graph measure distance features between the monozygotic twins, non-monozygotic siblings, and unrelated individuals. This includes Fig. 6c in the main text.
- ML_scripts/MEG_fMRI_GGM_LGM_ML_fingerprint.m (written in MATLAB 2020a) will run the genetic fingerprinting algorithm for classifying monozygotic twins, non-monozygotic siblings, and unrelated individuals.
- ML_scripts/plot_ML_fingerprint_results.m (written in MATLAB 2023b) will generate figures and tables of the genetic fingerprinting results. This includes Fig. 6b in the main text and Table S4-S7 in the Supplementary Information.

# Authors
Haatef Pourmotabbed
