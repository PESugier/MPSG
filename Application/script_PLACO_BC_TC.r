library(BhGLM)
library(arm)
library(tidyverse)
#library(Matrix)
library(tictoc)

library(MASS);library(arm);library(mvtnorm);library(invgamma);library(Matrix);
require(devtools)
source_url("https://github.com/RayDebashree/PLACO/blob/master/PLACO_v0.1.1.R?raw=TRUE")

# hard coded path directories
metapsg_dir <- "/home1/6_AMLAP/PE/METAPSG/scripts/" # path for mpsg scripts
datadir <- "/home1/6_AMLAP/PE/METAPSG/Applications/BC_TC_candidatepathways/data/" # data path
sumstatdir <- datadir # summary statistics path (generated from real data)
outputdir="/home1/6_AMLAP/PE/METAPSG/Applications/BC_TC_candidatepathways/app_like_gcpbayes_paper/results/metapsg/" # path for results
source(paste0(metapsg_dir,"metapsg_parallel.r"))

##############################################################################################################################
#
# Hardcoded parameters
#
no_cpu <- 60 # number of CPU for parallelisation
n_lambda_pts <- 30
verbose <- TRUE # Make it chatty

alpha_values <- seq(0.1,0.9,0.1)

##############################################################################################################################
# REAL Data (annotation info)
#
# load(paste0(datadir,"pleioData20191025_2300obs_3707snp.RData"))
# ls(data)
# Xtrue <- c()
# Ytrue <- c()
# Xb <- data$breast$Xb
# Xt <- data$thyroid$Xt
# # To Annot "GWAS" results
# genes <- unique(data$annot$gene)
# G <- ngrp <- ngenes <- length(genes)
# annoted_gene <- data$annot$gene
# annoted_snp <- as.data.frame(data$annot$snp)
# pathways <- unique(data$annot$pathway)
# npathways <- length(pathways)
# annoted_pathway <- data$annot$pathway
# # additional info
# n_s <- c(dim(Xb)[1],dim(Xt)[1])
# P <- dim(Xb)[2]
# save(genes, G, ngrp, ngenes, annoted_gene, annoted_snp, pathways, npathways, annoted_pathway, n_s, P, file=paste0(datadir,"info_annot_data_BC_TC.RData"))
load(paste0(datadir,"info_annot_data_BC_TC.RData"))


###########################################################
### TO PUT SUMMARY STATISTICS IN META-PSG OUTPUT FORMAT
load(paste0(sumstatdir,"summary_stat_application_gaussian_prior.RData"))
sumstat <-list()
###############################
Geneplotci <- Geneplotmed <- c()
#
for(g in 1:G){
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

###############################################
### RUNNING PLACO
##############################################################################################


  # Parameter PLACO
  pval_threshold <- 0.05

   sumstat <-list()
  ###############################
  #
  kt_g <- 0
  for(g in 1:G){
  #
    mg=length(Group[Group==g])
    sumstat[[g]] <- list(beta=NA, covar=NA, cov_inv=NA)
    sumstat[[g]]$beta <- array(NA,dim=c(mg,K))
    sumstat[[g]]$covar <- array(NA, dim=c(mg,mg,K))
    #
    for (k in 1:K){
    #
      sumstat[[g]]$beta[,k] <-  simBeta[k,Group==g]
      sumstat[[g]]$covar[,,k] <- SIGMA[[k]][Group==g,Group==g]
      ###################################
      snpnames <- 1:mg
      genename <- g
    }
  }

  Z1 <- sumstat[[1]]$beta[,1]/sqrt(diag(sumstat[[1]]$covar[,,1]))
  Z2 <- sumstat[[1]]$beta[,2]/sqrt(diag(sumstat[[1]]$covar[,,2]))
  for (g in 2:G){
    Z1 <- append(Z1,sumstat[[g]]$beta[,1]/sqrt(diag(sumstat[[g]]$covar[,,1])))
    Z2 <- append(Z2,sumstat[[g]]$beta[,2]/sqrt(diag(sumstat[[g]]$covar[,,2])))
  }
  pval1 <- (1-pnorm(abs(-Z1)))*2
  pval2 <- (1-pnorm(abs(-Z2)))*2
  # creation of Z and P matrices
  Z.matrix <- cbind(Z1, Z2)
  P.matrix <- cbind(pval1, pval2)
  k <- 2
  colnames(Z.matrix) <- paste("Z", 1:k, sep="")
  colnames(P.matrix) <- paste("P", 1:k, sep="")
  # performing a total correlation between two traits
  print(cor.test(Z1, Z2, method='pearson'))
  A <- cor.test(Z1, Z2, method='pearson')
  # deciding whether a "Decorrelate step" is needed or not  
  # based on comparing the obtained p_value between two traits and the "pval_threshold"
  if (A$p.value < pval_threshold) {
    writeLines("\n\n")
    print('It needs decorrelating the Z-scores')
    # Decorrelate the matrix of Z-scores
    R <- cor.pearson(Z.matrix, P.matrix, p.threshold=1e-4)
    # function for raising matrix to any power
    "%^%" <- function(x, pow)
      with(eigen(x), vectors %*% (values^pow * t(vectors)))
    Z.matrix.decor <- Z.matrix %*% (R %^% (-0.5))
    colnames(Z.matrix.decor) <- paste("ZD", 1:k, sep="")
  } else {
    writeLines("\n\n")
    print('No need to decorrelate the Z-scores')
    Z.matrix.decor <- Z.matrix
  }
  # ================================================================================
  # Obtaining the variance parameter estimates
  # ================================================================================  
  VarZ <- var.placo(Z.matrix.decor, P.matrix, p.threshold=1e-4)
  # ================================================================================
  # Applying test of pleiotropy for each variant
  # ================================================================================  
  out_1 <- sapply(1:P, function(i) placo(Z=Z.matrix.decor[i,], VarZ=VarZ))
  out_2 <- t(out_1)
  out_2 <- as.data.frame(out_2)
  out_2$T.placo <- as.numeric(out_2$T.placo)
  out_2$p.placo <- as.numeric(out_2$p.placo)
  # calculation of Adjusted p_value for each variant
  q.value <- p.adjust(out_2$p.placo, method='fdr', n=nrow(out_2))
  out_2 <- mutate(out_2, FDR.placo=q.value)
  ###############################################
  # Reporting method results for this simulation
  ##############################################################################################
  ResPleio_for_v[iseed,which(out_2$FDR.placo<=0.05)] <- 1 # Variables
  

