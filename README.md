# MPSG
Meta-analysis model adapted for Pleiotropy Selection with Group structure

Created by: Pierre-Emmanuel Sugier<br>
Creation date: 25 July 2024<br>
Update: July 2024<br>
<br>
<br>

**NOTES:**
<br><br>

This is a R package for the implementation of MPSG approach (Meta-analysis model adapted for Pleiotropy Selection with Group structure). 

This method performs penalized multivariate meta-analysis method adapted for pleiotropy and take into account the group structure information nested in the data to select relevant variants and genes (or pathways) from all the genetic information.


## Table of Contents
- [Method](#running-of-the-pipeline)
- [Application](#application)
- [Simulation study](#simulation-study)
- [Pipeline for pleiotropy analysis](#pipeline)
- [How to Cite](#how-to-cite)

## Method

The method performs penalized multivariate meta-analysis method adapted for pleiotropy and take into account the group structure information nested in the data to select relevant variants and genes (or pathways) by considering all genetic information simultaneously.

To do so, we implemented an alternating direction method of multipliers (ADMM) algorithm to perform both regularisation at variable level and group level.

The current release allows for:

1. Implementation of selection at group and variable level by using penalisation for multiple studies. If a variable or group of variables are selected, they are selected across all studies.

2. Choosing weights version for the penalisation of the method: standard version of the adaptive weights and no use of adaptive weights.

3. Parallelization along both hyperparameters alpha and lambda

4. Example code that includes illustrative bootstrapping code to detect variability in selected variables/groups of variables. 

TBD: (5. The use of the user friendly package for pleiotropy analysis, adapted to the use of MPSG from the pipeline of GCPBayes.)

Future versions may build out different versions of adaptive weights, and allow consideration of study selection in selected groups.


## Application

Our method is applied to the identification of potential pleiotropic genes between breast and thyroid cancer by using multivariate summary statistics generated from pathway candidate data.

## Simulation study

The performance of the novel approach is compared to benchmark gene-level and SNP-level meta-analysis approaches on simulated data by considering different kind of summary data as inputs.


## Pipeline for pleiotropy analysis 
TBD

## How to Cite
Sugier PE, Asgari Y, Sedki M, Truong T, Liquet B. "Meta-Analysis models with group structure for pleiotropy detection at gene and variant level by using summary statistics from multiple datasets"
<br>
<br>
