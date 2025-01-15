#!/usr/bin/env Rscript
library(optparse)

option_list = list(
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
require(devtools)
source_url("https://github.com/RayDebashree/PLACO/blob/master/PLACO_v0.1.1.R?raw=TRUE")

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

###############################################################################################################################

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

  snpnames <- 1:(G*mg)
  genename <- 1:G
  Geneplotci <- Geneplotmed <- c()

########################################################
# Running PLACO

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
  
  
  print(iseed)
}

Finaltime=Sys.time() - time00

############################################################################################################################################################################################
##############################################################################################
# Creating outputs
##############################################################################################

Summary_table_res_for_v <- data.frame(MCC=NA, MCCsd=NA, TPR=NA, TPRsd=NA, TNR=NA, TNRsd=NA)

Table_results_v <- data.frame(MCC=NA, TPR=NA, TNR=NA)


for(k in 1:nseed){
    Table_results_v[k,"TPR"] <- tfPositiveRate(ResPleio_for_v[k,], reals_v)$true_pos
    Table_results_v[k,"TNR"] <- 1-(tfPositiveRate(ResPleio_for_v[k,], reals_v)$false_pos)
    Table_results_v[k,"MCC"] <- tfPositiveRate(ResPleio_for_v[k,], reals_v)$MCC
}


# For variables
Summary_table_res_for_v$"MCC" <- mean(Table_results_v$"MCC")
Summary_table_res_for_v$"MCCsd" <- sd(Table_results_v$"MCC")
Summary_table_res_for_v$"TPR" <- mean(Table_results_v$"TPR")
Summary_table_res_for_v$"TPRsd" <- sd(Table_results_v$"TPR")
Summary_table_res_for_v$"TNR" <- mean(Table_results_v$"TNR")
Summary_table_res_for_v$"TNRsd" <- sd(Table_results_v$"TNR")


############################ 

res_summary_table <- Summary_table_res_for_v

# SAVE
write.table(res_summary_table, file=paste0(outputdir,"simu_res_placo.txt"), quote = FALSE)
write.table(Finaltime, file=paste0(outputdir,"time_placo.txt"), quote = FALSE)

