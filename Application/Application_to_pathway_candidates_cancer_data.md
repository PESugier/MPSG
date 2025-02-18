# Application to analysis of pleiotropy between thyroid cancer and breast cancer from pathway candidate summary data

In this vignette we used the following R packages

``` r
library(arm)
library(tidyverse)
library(doSnow)
library(tictoc)
```

For this example, we will use the R code version of the package. Hence, we need to source the main code that is sourcing the others. All codes need to be in the same folder.

``` r
source("metapsg_parallel.r")
```

## Hardcoded parameters

MPSG uses parallelisation to solve hyper-parameters tuning on a 2D-grid of alpha and lambda values. To do it, you need to specify the number of cpu you want to use for analysis. If you do not want to use parallellisation, you can directly use the metapsg_gv() function.

``` r
no_cpu <- 30 # number of CPU for parallelisation
```

The alpha parameter is controlling the variable/group penalisation ratio. As it is a value between 0 and 1, we suggest to test alpha for different 0.1 increasing values from 0.1 to 0.9.  The lambda parameter is controlling the overall penalisation. The user can control the number of lambda values tested. These values will be automatically generated following a logarithmic grid. We suggest to test for a 30 different lambda values per alpha value. 

``` r
alpha_values <- seq(0.1,0.9,0.1)
n_lambda_pts <- 30
```

## Load summary statistics and meta data

Then, we load precalculated summary statistics (here from a bayesian glm by using the bhglm packages), but you can for example use summary statistics from a GWAS. We also provide here SNP to genes and gene to pathways annotations for the data.

``` r
load("summary_stat_application_gaussian_prior.RData")
load("info_annot_data_BC_TC.RData")
```
