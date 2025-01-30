
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
head(annoted_snp)
```

Then, a correspondence vector (same size) which contains the name of the
gene (group) to which the SNP (variable/element) belongs

``` r
str(annoted_gene)
head(annoted_gene)
```

And finally, the vector of unique gene names

``` r
str(genes)
head(genes)
```

We also load from these data, the number of individuals for each study
from which the summary statistics has been generated, that is necessary
for the computation of the MIC.

``` r
n_s
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
sumstat <-list()
n_example <- 15
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

- here MPSG will compute adaptive weights as <code>adaptive=TRUE</code>

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

Pleiotropy at gene-level (group): two genes has been highlighted as
pleiotropic by MPSG

``` r
pleio_genes_names
```

Pleiotropy at SNP-level (variable): SNPs with non-zero values in
selected genes can be considered as pleiotropic

As example, the outputs of MPSG for the first gene are:

``` r
pleio_genes_names[1]
pleio_genes_res[[1]]
```

From this, we get the list of pleiotropic SNPs belonging to the gene.

``` r
names(which(pleio_genes_res[[1]][,1]!=0))
```
