# Simulation analysis
This section includes codes and files for the simulation analysis performed in the paper presenting the MPSG approach.

Created by: Pierre-Emmanuel Sugier<br>
Creation date: 25 July 2024<br>
Update: Jan 2025<br>
<br>

**NOTES:**
<br>

## Table of Contents
- [Simulation study](#simulation-study)
- [Codes](#codes)
- [Results](#results)

## Simulation study

The performance of the novel approach is compared to benchmark gene-level and SNP-level meta-analysis approaches on simulated data by considering different kind of summary data as inputs.

The approaches considered are:
- ASSET (SNP-level approach)
- PLACO (SNP-level approach)
- GCPBayes(DS function) (considering diagonal version and full version of covariance matrix)
- MPSG (considering diagonal version and full version of covariance matrix)
- MPSG (with adaptive weights) (considering diagonal version and full version of covariance matrix)

## Codes

Each simulation batch can be performed by running a Bash script with an initiation file informing on simulation parameters to consider, as follow: 

./00_Run_simulation_batch.sh simulation_batch_cor25_1.init

where simulation_batch_cor25_1.init is an example of .init file that has been run.

The Bash script is running several R scripts corresponding to each method considered one by one (ex: script_MetaPSG.R). The corresponding codes can be found in the Scripts/ folder.

## Results

Preprocessed results can be found in separated folders corresponding to different simulation batches.

Scripts used to generate these preprocessed results files from simulation outputs can also be found in this folder:
- script_generate_paper_tables_MCC.R
- script_generate_paper_tables_TNP.R
- script_generate_paper_tables_TPR.R
  
