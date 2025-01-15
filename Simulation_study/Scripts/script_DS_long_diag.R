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

pvalue=list()
ResPleio_for_g <- matrix(0,nseed,G)
ResPleio_for_v <- matrix(0,nseed,P)

ResBF <- ResFDR <- ResPleop <- matrix(0,nseed,G)
ResPleopSNPMEDIAN <- ResPleopSNPCI <- matrix(0,nseed,P)

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

# Running GCPBayes 
    for(g in 1:G){
        genename <- g
        snpnames <- seq(1,mg)
        list_betah <- list(simBeta[1,Group==g],simBeta[2,Group==g])
        list_sigmah <- list(diag(diag(SIGMA[[1]][Group==g,Group==g])),diag(diag(SIGMA[[2]][Group==g,Group==g])))
        # Launch of GCPBayes function (here DS)
        RES = DS(Betah = list_betah, Sigmah = list_sigmah, 
                             kappa0 = 0.5, 
                             sigma20 = 1,
                             m = mg,
                             K = K,
                             niter = 15000, 
                             burnin = 10000, 
                             nthin = 2,
                             nchains = 1, 
                             a1 = 0.1, 
                             a2 = 0.1, 
                             d1= 0.1, 
                             d2= 0.1,
                             snpnames = snpnames, 
                             genename = genename)
        Geneplotci <- append(Geneplotci,RES$Indicator$"Significant studies and Pleiotropic effect based on CI"$"Pleiotropic effect")
        Geneplotmed <- append(Geneplotmed,RES$Indicator$"Significant studies and Pleiotropic effect based on median thresholding"$"Pleiotropic effect")
        ResFDR[iseed,g] <- RES$Criteria$lBFDR
        ResPleop[iseed,g] <- RES$Criteria$theta
        ResBF[iseed,g] <- RES$Criteria$log10BF
    }
  ResPleopSNPCI[iseed,]=Geneplotci
  ResPleopSNPMEDIAN[iseed,]=Geneplotmed
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

model_g <- matrix(0,nseed,G)
#trbf=c(1,2);
trf=0.05
theta=0.5
for(k in 1:nseed){
  for(j in 1:G){
    #if((ResBF[k,j]>trbf[1]) & (ResPleop[k,j]>0.5))(modelb1[k,j]=1)
    #if((ResBF[k,j]>trbf[2]) & (ResPleop[k,j]>0.5))(modelb2[k,j]=1)
    if((ResFDR[k,j]<trf) & (ResPleop[k,j]>theta))(model_g[k,j] <- 1)
  }}

for(k in 1:nseed){
  Table_results_g[k,"TPR"] <- tfPositiveRate(model_g[k,], reals_g)$true_pos
  Table_results_g[k,"TNR"] <- 1 - tfPositiveRate(model_g[k,], reals_g)$false_pos
  Table_results_g[k,"MCC"] <- tfPositiveRate(model_g[k,], reals_g)$MCC
}

################################## SNP- Pleiotropy##################
RESPlEOPsnp=RRgene=matrix(0,nseed,P)
for(k in 1:nseed){
  M1=rep(model_g[k,],each=mg)
  for(j in 1:P){
    if((ResPleopSNPCI[k,j]=="Yes") & (M1[j]==1))(RESPlEOPsnp[k,j]=1)
  }}

############################ 
# reals=rep(0,P)
# reals[BETA2!=0]=1

for(k in 1:nseed){
    Table_results_v[k,"TPR"] <- tfPositiveRate(RESPlEOPsnp[k,], reals_v)$true_pos
    Table_results_v[k,"TNR"] <- 1-(tfPositiveRate(RESPlEOPsnp[k,], reals_v)$false_pos)
    Table_results_v[k,"MCC"] <- tfPositiveRate(RESPlEOPsnp[k,], reals_v)$MCC
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
write.table(res_summary_table, file=paste0(outputdir,"simu_res_ds_long_diag.txt"), quote = FALSE)
write.table(Finaltime, file=paste0(outputdir,"time_ds_long_diag.txt"), quote = FALSE)

