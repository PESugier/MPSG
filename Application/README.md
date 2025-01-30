
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Example of using MPSG approach for analysis of pleiotropy: by using a gene selection from thyroid cancer and breast cancer summary data

<!-- badges: start -->
<!-- badges: end -->

In this vignette we used the following R packages

``` r
library(Matrix)
library(Rcpp)
library(RcppArmadillo)
library(arm)
library(tictoc)
library(parallel)
library(doParallel)
library(doSNOW)
library(abind)
```

For this example, we will use the R code version of the package. Hence,
we need to source the main code that is sourcing the others. All codes
need to be in the same folder.

``` r
source("metapsg_parallel.r")
```

## Hardcoded parameters

MPSG uses parallelisation to solve hyper-parameters tuning on a 2D-grid
of alpha and lambda values. To do it, you need to specify the number of
cpu you want to use for analysis. If you do not want to use
parallellisation, you can directly use the metapsg_gv() function. For a
fast running of the method, we recommend to set up the number of cpu
equal to the number of values of the lambda grid (or half of it).

``` r
no_cpu <- 5 # number of CPU for parallelisation
```

The alpha parameter is controlling the variable/group penalisation
ratio. As it is a value between 0 and 1, we suggest to test for alpha
values between 0.1 to 0.9 in steps of 0.1. For this example, we will use
0.2 steps. The lambda parameter is controlling the overall penalisation.
The user can control the number of lambda values tested. These values
will be automatically generated following a logarithmic grid. We suggest
to test for a 30 different lambda values per alpha value. For this
example, we will use 10 different lambda values.

``` r
alpha_values <- seq(0.1,0.9,0.2)
n_lambda_pts <- 10
```

## Load summary statistics and meta data

We load precalculated summary statistics (here from a bayesian glm by
using the bhglm packages), but you can for example use summary
statistics from a GWAS. We also provide here SNPs-to-genes and
gene-to-pathways annotations for the corresponding data, by loading the
info_annot_data_BC_TC.RData file.

``` r
# load summary statistics
load("summary_stat_application_gaussian_prior.RData")
# load annotation file
load("info_annot_data_BC_TC.RData")
```

First, a single column data.frame that contain names of each SNP
(variable/element)

``` r
str(annoted_snp)
#> 'data.frame':    3707 obs. of  1 variable:
#>  $ data$annot$snp: chr  "rs7531583" "rs12044597" "rs697679" "rs697680" ...
```

``` r
head(annoted_snp)
#>   data$annot$snp
#> 1      rs7531583
#> 2     rs12044597
#> 3       rs697679
#> 4       rs697680
#> 5       rs228729
#> 6       rs875994
```

Then, a correspondence vector (same size) which contains the name of the
gene (group) to which the SNP (variable/element) belongs

``` r
str(annoted_gene)
#>  chr [1:3707] "NADK" "NADK" "PER3" "PER3" "PER3" "PER3" "PER3" "PER3" "PER3" "PER3" "PER3" "PER3" "PER3" "PER3" "NMNAT1" "NMNAT1" ...
```

``` r
head(annoted_gene)
#> [1] "NADK" "NADK" "PER3" "PER3" "PER3" "PER3"
```

And finally, the vector of unique gene names

``` r
str(genes)
#>  chr [1:331] "NADK" "PER3" "NMNAT1" "MTHFR" "FUCA1" "NR0B2" "RPA2" "NT5C1A" "MUTYH" "AKR1A1" "LEPR" "NEGR1" "TGFBR3" "GSTM4" ...
```

``` r
head(genes)
#> [1] "NADK"   "PER3"   "NMNAT1" "MTHFR"  "FUCA1"  "NR0B2"
```

We also load from these data, the number of individuals for each study
from which the summary statistics has been generated, that is necessary
for the computation of the MIC.

``` r
n_s
#> [1] 2297 2303
```

## Preparation of MPSG inputs

The MPSG approach require a list where each element is corresponding to
a group, here the genes (that could be also pathways for example). In
each group: -
<code>$beta</code> will contain a matrix of beta estimates: each column corresponding to a separate study, and each line to a variable, element of the group.
- <code>$covar</code> will contain covariance matrices for each study:
<code>$covar[,,1]</code> for study 1, <code>covar[,,2]</code> for study 2, etc.
- then, <code>$cov_inv</code> will contain matrices that can be
calculated from \$covar by using the <code>f_cpp_inv_var()</code>
function of the MPSG package.

For example purpose, we will only use the first 15 genes of the panel
for analysis.

``` r
n_example <- 15
genes[1:n_example]
#>  [1] "NADK"   "PER3"   "NMNAT1" "MTHFR"  "FUCA1"  "NR0B2"  "RPA2"   "NT5C1A" "MUTYH"  "AKR1A1" "LEPR"   "NEGR1"  "TGFBR3" "GSTM4" 
#> [15] "GSTM2"
```

``` r
sumstat <-list()
#
for(g in 1:n_example){ # for example purpose
#
  # print(g)
  selected_ids = which(annoted_gene==genes[g])
  mg <- length(selected_ids)
  sumstat[[g]] <- list(beta=NA, covar=NA, cov_inv=NA)
  sumstat[[g]]$beta <- array(NA,dim=c(mg,2))
  sumstat[[g]]$covar <- array(0, dim=c(mg,mg,2))
  sumstat[[g]]$cov_inv <- array(NA, dim=c(mg,mg,2))
  #
  sumstat[[g]]$beta[,1] <-  SUMSTAT_bglm_gaussian[[g]]$betah1
  sumstat[[g]]$beta[,2] <-  SUMSTAT_bglm_gaussian[[g]]$betah2
  rownames(sumstat[[g]]$beta) <- annoted_snp[selected_ids,]
  sumstat[[g]]$covar[,,1] <- SUMSTAT_bglm_gaussian[[g]]$sigmah1
  sumstat[[g]]$covar[,,2] <- SUMSTAT_bglm_gaussian[[g]]$sigmah2
  tryCatch({
    sumstat[[g]]$cov_inv <- f_cpp_inv_var(covar=sumstat[[g]]$covar, same_var=FALSE)
  },error = function(e){
    print(paste0("Une erreur est survenue pour le gène numéro: ", g," -> ", e$message))
    sumstat[[g]]$covar[,,1] <<- sumstat[[g]]$covar[,,1] + 0.01*diag(diag(sumstat[[g]]$covar[,,1]))
    sumstat[[g]]$covar[,,2] <<- sumstat[[g]]$covar[,,2] + 0.01*diag(diag(sumstat[[g]]$covar[,,2]))
    sumstat[[g]]$cov_inv <<- f_cpp_inv_var(covar=sumstat[[g]]$covar, same_var=FALSE)
    print("Erreur corrigée avec l'ajout d'un petit terme diagonal")
  })
}
```

## Run MPSG

We can now run the MPSG method on <code>sumstat</code> inputs, by using
the <code>f_selec_mod_metapsg_gv_parallel()</code> function that.

- here MPSG will be run with adaptive weights as
  <code>adaptive=TRUE</code>

``` r
res <- f_selec_mod_metapsg_gv_parallel(sumstat, n_s, alpha=alpha_values, adaptive=TRUE, same_var = FALSE, rho=1, t_inc=2, t_dec=2, mu=10, max_iter = 1000, tol_abs = 1e-03, tol_rel = 1e-02, eps=1e-03, npts=n_lambda_pts, ncpu=no_cpu)
```

Then, we annotate and save the selected model

``` r
final_model <- res$model
names(final_model) <- genes[1:n_example]
for (g in 1:n_example){
  selected_ids = which(annoted_gene==genes[g])
  rownames(final_model[[g]]) <- annoted_snp[selected_ids,]
}
```

We select names and indices of genes with non-zero elements.

``` r
pleio_genes_res <- list()
pleio_genes_names <- c()
pleio_genes_indice <- c()
for (g in 1:n_example){
  if(norm(final_model[[g]])!=0){
    pleio_genes_res <- append(pleio_genes,list(c=final_model[[g]]))
    pleio_genes_names <- c(pleio_genes_names,names(final_model)[g])
    pleio_genes_indice <- c(pleio_genes_indice,g)
  }
}
```

## Interpretation of the results

### Pleiotropy at gene-level (group)

Two genes has been highlighted as pleiotropic by MPSG

``` r
pleio_genes_names
#> [1] "NEGR1"  "TGFBR3"
```

### Pleiotropy at SNP-level (variable)

SNPs with non-zero values in selected genes can be considered as
pleiotropic

As example, the outputs of MPSG for the first gene are:

``` r
pleio_genes_names[1]
#> [1] "NEGR1"
```

``` r
pleio_genes_res[[1]]
#>                    [,1]         [,2]
#> rs1506457   0.000000000  0.000000000
#> rs12142099 -0.035124321 -0.021367845
#> rs17091209  0.000000000  0.000000000
#> rs7519228   0.000000000  0.000000000
#> rs12126700  0.000000000  0.000000000
#> rs7544223   0.000000000  0.000000000
#> rs12135145  0.000000000  0.000000000
#> rs2056206   0.000000000  0.000000000
#> rs12142367  0.000000000  0.000000000
#> rs1517756   0.000000000  0.000000000
#> rs12142140  0.000000000  0.000000000
#> rs357221    0.000000000  0.000000000
#> rs1889298   0.000000000  0.000000000
#> rs12129849 -0.055514003 -0.158768014
#> rs17503667  0.000000000  0.000000000
#> rs17569138  0.000000000  0.000000000
#> rs17503709  0.064511262 -0.027287791
#> rs357233    0.000000000  0.000000000
#> rs357214   -0.079548491 -0.049896615
#> rs11209800  0.110212725 -0.004108690
#> rs928615    0.000000000  0.000000000
#> rs1762732   0.000000000  0.000000000
#> rs6696052   0.000000000  0.000000000
#> rs11209806 -0.139833542  0.048346762
#> rs12143094  0.000000000  0.000000000
#> rs11578555  0.196881375 -0.151611128
#> rs7524126   0.000000000  0.000000000
#> rs591031   -0.071292695  0.143199539
#> rs2422014   0.000000000  0.000000000
#> rs11209815  0.000000000  0.000000000
#> rs11209817  0.000000000  0.000000000
#> rs1591382   0.000000000  0.000000000
#> rs10889930  0.000000000  0.000000000
#> rs12023733  0.000000000  0.000000000
#> rs10789324 -0.219180832 -0.112692209
#> rs2125777  -0.001846089 -0.011752121
#> rs4649946   0.000000000  0.000000000
#> rs17091734  0.049827961 -0.045149620
#> rs1459802   0.000000000  0.000000000
#> rs531071    0.000000000  0.000000000
#> rs12031844  0.086216746  0.343746182
#> rs12141391  0.000000000  0.000000000
#> rs12144209 -0.082546902  0.178950608
#> rs6424449   0.000000000  0.000000000
#> rs1318907   0.010093536 -0.094355557
#> rs12119234  0.000000000  0.000000000
#> rs1426178   0.000000000  0.000000000
#> rs1426177   0.000000000  0.000000000
#> rs574888    0.000000000  0.000000000
#> rs11579239  0.000000000  0.000000000
#> rs518523    0.000000000  0.000000000
#> rs17091876  0.000000000  0.000000000
#> rs601452    0.000000000  0.000000000
#> rs2114214   0.000000000  0.000000000
#> rs12133406  0.000000000  0.000000000
#> rs2114213   0.000000000  0.000000000
#> rs17091962  0.000000000  0.000000000
#> rs6703903   0.000000000  0.000000000
#> rs1937188   0.000000000  0.000000000
#> rs12133119 -0.312613982  0.163547869
#> rs6696193  -0.052069882 -0.042196475
#> rs12566464  0.000000000  0.000000000
#> rs6661979   0.000000000  0.000000000
#> rs12139692  0.000000000  0.000000000
#> rs1954600   0.141877286 -0.041496058
#> rs1872578   0.000000000  0.000000000
#> rs17363292  0.000000000  0.000000000
#> rs1032154   0.000000000  0.000000000
#> rs12119337  0.000000000  0.000000000
#> rs11806908 -0.208748025 -0.097182761
#> rs1486090   0.000000000  0.000000000
#> rs1155261   0.000000000  0.000000000
#> rs1486087   0.000000000  0.000000000
#> rs2821278   0.000000000  0.000000000
#> rs1486086   0.000000000  0.000000000
#> rs4348675   0.000000000  0.000000000
#> rs2200485   0.000000000  0.000000000
#> rs12409966  0.000000000  0.000000000
#> rs1486105  -0.207300886  0.068261607
#> rs12091740  0.000000000  0.000000000
#> rs2821248   0.000000000  0.000000000
#> rs12128707  0.000000000  0.000000000
#> rs17588812  0.033802059  0.124878454
#> rs11587434 -0.302031134  0.087485186
#> rs1157072   0.000000000  0.000000000
#> rs1342862  -0.009067584 -0.031930285
#> rs6656980   0.000000000  0.000000000
#> rs2821285  -0.075250544  0.058405226
#> rs1776012   0.000000000  0.000000000
#> rs12127789  0.000000000  0.000000000
#> rs9326098  -0.016685054 -0.009650582
#> rs3101338   0.000000000  0.000000000
#> rs3101336   0.000000000  0.000000000
```

From this, we get the list of pleiotropic SNPs belonging to the gene.

``` r
names(which(pleio_genes_res[[1]][,1]!=0))
#>  [1] "rs12142099" "rs12129849" "rs17503709" "rs357214"   "rs11209800" "rs11209806" "rs11578555" "rs591031"   "rs10789324" "rs2125777" 
#> [11] "rs17091734" "rs12031844" "rs12144209" "rs1318907"  "rs12133119" "rs6696193"  "rs1954600"  "rs11806908" "rs1486105"  "rs17588812"
#> [21] "rs11587434" "rs1342862"  "rs2821285"  "rs9326098"
```
