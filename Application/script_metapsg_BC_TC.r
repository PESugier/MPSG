library(arm)
library(tidyverse)
#library(Matrix)
library(tictoc)

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


##############################################
## METAPSG GV SOLVER WITH ADAPTIVE WEIGHTS
#############################################################################################
n_alphas <- length(alpha_values)
t1 <- Sys.time()
res <- f_selec_mod_metapsg_gv_parallel(sumstat, n_s, alpha=alpha_values, adaptive=TRUE, same_var = FALSE, rho=1, t_inc=2, t_dec=2, mu=10, max_iter = 1000, tol_abs = 1e-03, tol_rel = 1e-02, eps=1e-03, npts=n_lambda_pts, ncpu=no_cpu)
# save MIC values
t2 <- Sys.time()
print(t2 - t1)
########## SAVE MODEL ############
final_model <- res$model
names(final_model) <- genes
for (g in 1:G){
  selected_ids = which(annoted_gene==genes[g])
  rownames(final_model[[g]]) <- annoted_snp[selected_ids,]
}
##############################################################################################
RESULT FILE MANAGEMENT
################################
pleio_genes<-list()
indice_hit<-c()
for (g in 1:G){
  if(norm(final_model[[g]])!=0){
    pleio_genes <- append(pleio_genes,list(c=final_model[[g]]))
    indice_hit <- c(indice_hit,g)
  }
}
names(pleio_genes) <- genes[indice_hit]
metapsg_pleiotropic_genes <- pleio_genes
save(metapsg_pleiotropic_genes, file=paste0(outputdir, "RESULTS_MetaPSG_Best_Model_adaptive_weights.RData"))


###############################################
### METAPSG GV SOLVER WITHOUT ADAPTIVE WEIGHTS
##############################################################################################
n_alphas <- length(alpha_values)
t1 <- Sys.time()
res <- f_selec_mod_metapsg_gv_parallel(sumstat, n_s, alpha=alpha_values, adaptive=FALSE, same_var = FALSE, rho=1, t_inc=2, t_dec=2, mu=10, max_iter = 1000, tol_abs = 1e-03, tol_rel = 1e-02, eps=1e-03, npts=n_lambda_pts, ncpu=no_cpu)
# save MIC values
t2 <- Sys.time()
print(t2 - t1)
########## SAVE MODEL ############
final_model <- res$model
names(final_model) <- genes
for (g in 1:G){
  selected_ids = which(annoted_gene==genes[g])
  rownames(final_model[[g]]) <- annoted_snp[selected_ids,]
}
###############################################################################################
# RESULT FILE MANAGEMENT
#################################
pleio_genes<-list()
indice_hit<-c()
for (g in 1:G){
  if(norm(final_model[[g]])!=0){
    pleio_genes <- append(pleio_genes,list(c=final_model[[g]]))
    indice_hit <- c(indice_hit,g)
  }
}
names(pleio_genes) <- genes[indice_hit]
metapsg_pleiotropic_genes <- pleio_genes
save(metapsg_pleiotropic_genes, file=paste0(outputdir, "RESULTS_MetaPSG_Best_Model.RData"))

