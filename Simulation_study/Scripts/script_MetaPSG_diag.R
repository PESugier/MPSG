#!/usr/bin/env Rscript
library(optparse)

option_list = list(
  make_option(c("--path_metapsg"), type="character", default=NULL, 
              help="working directory", metavar="character"),
  make_option(c("--path_out"), type="character", default=NULL, 
              help="generic input file names", metavar="character"),
  make_option(c("--n_rep"), type="numeric", default=NULL, 
              help="number of replications", metavar="number"),
  make_option(c("--effect"), type="numeric", default=NULL, 
              help="Amount of variable effect", metavar="number"),
  make_option(c("--SNR"), type="numeric", default=NULL, 
              help="signal to noize ratio", metavar="number"),
  make_option(c("--GS"), type="numeric", default=NULL, 
              help="percentage of group sparsity", metavar="number"),
  make_option(c("--IGS"), type="numeric", default=NULL, 
              help="percentage of sparsity for variable in groups", metavar="number"),
  make_option(c("--n_study"), type="numeric", default=NULL, 
              help="number of studies considered", metavar="number"),
  make_option(c("--n_indiv"), type="numeric", default=NULL, 
              help="number of simulated individuals in each study", metavar="number"),
  make_option(c("--g"), type="numeric", default=NULL, 
              help="number of simulated groups per replication", metavar="numeric"),
  make_option(c("--mg"), type="numeric", default=NULL, 
              help="number of variables per group", metavar="numeric"),
  make_option(c("--corr"), type="numeric", default=NULL, 
              help="correlation between variables in the same group", metavar="numeric")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# Load other required libraries
library(MASS);library(arm);library(mvtnorm);library(invgamma);library(Matrix);library(GCPBayes)
library(BhGLM)
library(arm)
library(dplyr)

# load Meta-PSG functions
metadir <- opt$path_metapsg
source(paste0(metadir,"metapsg_gv_solver.r"))
source(paste0(metadir,"metapsg_gsv_solver.r"))
source(paste0(metadir,"metapsg_parallel.r"))

time00 = Sys.time()


##################################################
# PATH
outputdir <- opt$path_out

##################################################
# PARAMETERS

nseed <- opt$n_rep
K <- opt$n_study # Number of studies
mg <- opt$mg
G <- ngrp <- opt$g
P <- mg*G
var_corr <- opt$corr

# Signal / Noize
effect <- opt$effect
SNR <- opt$SNR

# Sparsity intra/inter groups
IGS <- opt$IGS
intraGroupSparsity <- IGS/100
groupSparsity <- opt$GS/100

# N: Effectives
N <- opt$n_indiv # A VERIFIER
n_s <- rep(N,K)

reals_g <- rep(0,G)
reals_g <- c(rep(1,G*(1-groupSparsity)),rep(0,G*groupSparsity))

reals_v <- rep(0,G*mg)
structure_v_in_true_g <- c(rep(1, mg*(1-intraGroupSparsity)),rep(0,mg*intraGroupSparsity))
structure_v_in_false_g <- rep(0,mg)
reals_v <- c(rep(structure_v_in_true_g,G*(1-groupSparsity)),rep(structure_v_in_false_g,G*groupSparsity))

###############################################################################################################################

# helping function 
tfPositiveRate = function(est_model, true_model) {
    n_true_pos = sum(est_model * true_model)
    n_true_model = sum(true_model)
    n_par = length(c(true_model))
    n_pos = sum(est_model)
    # True/False positive rate
    TP <- n_true_pos
    FP <- n_pos - n_true_pos
    true_pos = n_true_pos / n_true_model
    false_pos = (n_pos - n_true_pos) / (n_par - n_true_model)
    FN <- n_true_model - n_true_pos
    TN <- length(true_model) - n_true_model - (n_pos - n_true_pos)
    N <- TN+TP+FN+FP
    P=(TP+FP)/N
    S=(TP+FN)/N
    test <- sqrt(P*S*(1-S)*(1-P))
    if (test==0) MCC <- 0 else MCC <- (TP/N-S*P)/sqrt(P*S*(1-S)*(1-P))
    return(list(true_pos = true_pos, false_pos = false_pos,MCC=MCC))
}

###############################################################################################################################
# To generate the SIGMA covariance matrix of beta estimates from simulated X data
# (generated for one group & reproduced once for all group of both studies)

group = parallel::splitIndices(P, ngrp) # groups as a list of indices
Group=as.numeric(gl(G, P/G))

# simulating correlated markers values through an additive genetic model
# X <- scale(as.matrix(simX), center = TRUE, scale = FALSE)

###############################################################################################################################

pvalue=list()
ResPleio_for_g <- matrix(0,nseed,G)
ResPleio_for_v <- matrix(0,nseed,P)


# Simulation loop
for (iseed in 1:nseed){
  set.seed(iseed)
  
  # Simulation of cov(X)
  R <- matrix(var_corr,mg,mg)+(1-var_corr)*diag(mg)
  CovX <- kronecker(diag(G), R)
  
  # Betas
  Group=as.numeric(gl(G, mg))
  sign=rbinom(mg,1,.5)
  sign[sign==0]=-1
  betas=rep(0,mg)
  betat <- c(rep(effect,mg*(1-intraGroupSparsity)),rep(0,mg*intraGroupSparsity))
  
  # True Betas
  betat=betat*sign
  BETA1=c(rep(betat,(1-groupSparsity)*G),rep(betas,groupSparsity*G))
  BETA2=BETA1
  BETA2[Group==2]=-BETA2[Group==2]
  
  # SIGMA
  noize1 <- t(BETA1)%*% CovX %*% BETA1 / SNR
  noize2 <- t(BETA2)%*% CovX %*% BETA2 / SNR
  SIGMA <- list()
  SIGMA[[1]] <- as.numeric(noize1/N)*solve(CovX)
  SIGMA[[2]] <- as.numeric(noize2/N)*solve(CovX)
  
  simBeta <- rbind(rmvnorm(1,BETA1,SIGMA[[1]]),rmvnorm(1,BETA2,SIGMA[[2]]))
  #SIGMA=diag(diag(SIGMA))

  snpnames=1:(G*mg)
  genename=1:G

#################################################################
# Inputs for MetaPSG

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
      sumstat[[g]]$cov_inv <- array(NA, dim=c(mg,mg,K))
      #
      for (k in 1:K){
      #
        sumstat[[g]]$beta[,k] <-  simBeta[k,Group==g]
        sumstat[[g]]$covar[,,k] <- diag(diag(SIGMA[[k]][Group==g,Group==g]))
        ###################################
        snpnames <- 1:mg
        genename <- g
      }
      sumstat[[g]]$cov_inv <- f_cpp_inv_var(covar=sumstat[[g]]$covar, same_var=FALSE)
    }

#################################################################
# Running MetaPSG
    Res_grid_alpha <- list()
    alpha_values <- seq(0.1,0.9,0.1)
    n_alphas <- length(alpha_values)
    n_lambda_pts <- 50
    MIC_grid <- matrix(NA, nrow=n_lambda_pts, ncol=n_alphas)
    t1 <- Sys.time()
    for (i in 1:n_alphas){
      # To generate results for this value of alpha
      Res_grid_alpha[[i]] <- f_selec_mod_metapsg_gv(sumstat, n_s, alpha=alpha_values[i], adaptive=FALSE, t_inc=2, npts=n_lambda_pts)
      # save MIC values
      MIC_grid[,i] <- Res_grid_alpha[[i]]$selmod$MIC
    }
    t2 <- Sys.time()
    print(t2 - t1)
    # index of the Best Model according to MIC:
    min_index <- which(MIC_grid == min(MIC_grid), arr.ind = TRUE)
    if(length(min_index)!=2){min_index <- min_index[1,]} # In case of more than one model minimizing the MIC (very likely to be both the same, with all estimates = 0...)
    # Best Model:
    BestModel <- Res_grid_alpha[[min_index[2]]]$results[[min_index[1]]]$beta_select

    
    ###############################################
    # Reporting method results for this simulation
    ##############################################################################################
    ResPleio_for_g[iseed,which(sapply(FUN=norm, BestModel, simplify=TRUE)!=0)] <- 1 # Groups
    #
    BM_v <- BestModel
    for (i in 1:G){
      BM_v[[i]] <- BM_v[[i]]^2
      BM_v[[i]] <- apply(BM_v[[i]], MARGIN = 1, FUN = sum) %>% sqrt
    }
    BM_v <- unlist(BM_v)
    ResPleio_for_v[iseed,which(BM_v!=0)] <- 1 # Variables

    print(iseed)
}

Finaltime=Sys.time() - time00

############################################################################################################################################################################################
##############################################################################################
# Creating outputs
##############################################################################################

Summary_table_res_for_g <- data.frame(MCC=NA, MCCsd=NA, TPR=NA, TPRsd=NA, TNR=NA, TNRsd=NA)
Summary_table_res_for_v <- data.frame(MCC=NA, MCCsd=NA, TPR=NA, TPRsd=NA, TNR=NA, TNRsd=NA)

Table_results_g <- data.frame(MCC=NA, TPR=NA, TNR=NA)
Table_results_v <- data.frame(MCC=NA, TPR=NA, TNR=NA)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Pleiotropy-total @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  


for(k in 1:nseed){
  Table_results_g[k,"TPR"] <- tfPositiveRate(ResPleio_for_g[k,], reals_g)$true_pos
  Table_results_g[k,"TNR"] <- 1 - tfPositiveRate(ResPleio_for_g[k,], reals_g)$false_pos
  Table_results_g[k,"MCC"] <- tfPositiveRate(ResPleio_for_g[k,], reals_g)$MCC
}

for(k in 1:nseed){
    Table_results_v[k,"TPR"] <- tfPositiveRate(ResPleio_for_v[k,], reals_v)$true_pos
    Table_results_v[k,"TNR"] <- 1-(tfPositiveRate(ResPleio_for_v[k,], reals_v)$false_pos)
    Table_results_v[k,"MCC"] <- tfPositiveRate(ResPleio_for_v[k,], reals_v)$MCC
}


# For groups
Summary_table_res_for_g$"MCC" <- mean(Table_results_g$"MCC")
Summary_table_res_for_g$"MCCsd" <- sd(Table_results_g$"MCC")
Summary_table_res_for_g$"TPR" <- mean(Table_results_g$"TPR")
Summary_table_res_for_g$"TPRsd" <- sd(Table_results_g$"TPR")
Summary_table_res_for_g$"TNR" <- mean(Table_results_g$"TNR")
Summary_table_res_for_g$"TNRsd" <- sd(Table_results_g$"TNR")
# For variables
Summary_table_res_for_v$"MCC" <- mean(Table_results_v$"MCC")
Summary_table_res_for_v$"MCCsd" <- sd(Table_results_v$"MCC")
Summary_table_res_for_v$"TPR" <- mean(Table_results_v$"TPR")
Summary_table_res_for_v$"TPRsd" <- sd(Table_results_v$"TPR")
Summary_table_res_for_v$"TNR" <- mean(Table_results_v$"TNR")
Summary_table_res_for_v$"TNRsd" <- sd(Table_results_v$"TNR")


############################ 

rownames(Summary_table_res_for_g) <- paste0("g_", effect)
rownames(Summary_table_res_for_v) <- paste0("v_", effect)
res_summary_table <- rbind(Summary_table_res_for_g, Summary_table_res_for_v)

# SAVE
write.table(res_summary_table, file=paste0(outputdir,"simu_res_metapsg_diag.txt"), quote = FALSE)
write.table(Finaltime, file=paste0(outputdir,"time_metapsg_diag.txt"), quote = FALSE)

