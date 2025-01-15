# functions

metapsgdir="/home1/6_AMLAP/PE/METAPSG/scripts/"
scriptsimudir="/home1/6_AMLAP/PE/METAPSG/SIMU/beta_from_CovX/"

f_create_table_groupe <- function(names_folders, name_out){

  v_N <- c("1k"," ","1.5k"," ","2k"," ","2.5k"," ","3k", " ")
  v_D_ND <- c("ND","D","ND","D","ND","D","ND","D","ND","D")
  MCC_g <- data.frame(N=v_N, inputs=v_D_ND, DS=NA, DS_v=NA, PSG=NA, PSG_v=NA, PSG_aw=NA, PSG_aw_v=NA)

  for (i in 1:length(names_folders)){
    name_batch=names_folders[i]
    outputdir=paste0(resdir,name_batch)
    placo=read.table(paste0(outputdir,"simu_res_placo.txt"))
    asset=read.table(paste0(outputdir,"simu_res_asset.txt"))
    metapsg=read.table(paste0(outputdir,"simu_res_metapsg.txt"))
    metapsg_diag=read.table(paste0(outputdir,"simu_res_metapsg_diag.txt"))
    metapsg_adapt=read.table(paste0(outputdir,"simu_res_metapsg_adapt.txt"))
    metapsg_diag_adapt=read.table(paste0(outputdir,"simu_res_metapsg_diag_adapt.txt"))
    ds=read.table(paste0(outputdir,"simu_res_ds_short.txt"))
    ds_diag=read.table(paste0(outputdir,"simu_res_ds_short_diag.txt"))
    # GROUPES
    #
    MCC_g$DS[which(MCC_g$inputs=="ND")][i] <- round(ds$MCC[1], digits=2)
    MCC_g$DS[which(MCC_g$inputs=="D")][i] <- round(ds_diag$MCC[1], digits=2)
    #    
    MCC_g$DS_v[which(MCC_g$inputs=="ND")][i] <- round(ds$MCCsd[1], digits=2)
    MCC_g$DS_v[which(MCC_g$inputs=="D")][i] <- round(ds_diag$MCCsd[1], digits=2)
    #    
    MCC_g$PSG[which(MCC_g$inputs=="ND")][i] <- round(metapsg$MCC[1], digits=2)
    MCC_g$PSG[which(MCC_g$inputs=="D")][i] <- round(metapsg_diag$MCC[1], digits=2)
    #    
    MCC_g$PSG_v[which(MCC_g$inputs=="ND")][i] <- round(metapsg$MCCsd[1], digits=2)
    MCC_g$PSG_v[which(MCC_g$inputs=="D")][i] <- round(metapsg_diag$MCCsd[1], digits=2)
    #    
    MCC_g$PSG_aw[which(MCC_g$inputs=="ND")][i] <- round(metapsg_adapt$MCC[1], digits=2)
    MCC_g$PSG_aw[which(MCC_g$inputs=="D")][i] <- round(metapsg_diag_adapt$MCC[1], digits=2)
    #
    MCC_g$PSG_aw_v[which(MCC_g$inputs=="ND")][i] <- round(metapsg_adapt$MCCsd[1], digits=2)
    MCC_g$PSG_aw_v[which(MCC_g$inputs=="D")][i] <- round(metapsg_diag_adapt$MCCsd[1], digits=2)
  }
  write.table(MCC_g, file=paste0(resdir,name_out,"_groupe.txt"), quote = FALSE, sep = " &\ ")
}

f_create_table_variable <- function(names_folders, name_out){
  #
  v_N <- c("1k"," ","1.5k"," ","2k"," ","2.5k"," ","3k", " ")
  v_D_ND <- c("ND","D","ND","D","ND","D","ND","D","ND","D")
  MCC_v <- data.frame(N=v_N, inputs=v_D_ND, DS=NA, DS_v=NA, PSG=NA, PSG_v=NA, PSG_aw=NA, PSG_aw_v=NA, ASSET=NA, ASSET_v=NA, PLACO=NA, PLACO_v=NA)
  #
  for (i in 1:length(names_folders)){
    name_batch=names_folders[i]
	outputdir=paste0(resdir,name_batch)
    placo=read.table(paste0(outputdir,"simu_res_placo.txt"))
    asset=read.table(paste0(outputdir,"simu_res_asset.txt"))
    metapsg=read.table(paste0(outputdir,"simu_res_metapsg.txt"))
    metapsg_diag=read.table(paste0(outputdir,"simu_res_metapsg_diag.txt"))
    metapsg_adapt=read.table(paste0(outputdir,"simu_res_metapsg_adapt.txt"))
    metapsg_diag_adapt=read.table(paste0(outputdir,"simu_res_metapsg_diag_adapt.txt"))
    ds=read.table(paste0(outputdir,"simu_res_ds_short.txt"))
    ds_diag=read.table(paste0(outputdir,"simu_res_ds_short_diag.txt"))
    # variables
	#
    MCC_v$DS[which(MCC_v$inputs=="ND")][i] <- round(ds$MCC[2], digits=2)
    MCC_v$DS[which(MCC_v$inputs=="D")][i] <- round(ds_diag$MCC[2], digits=2)
    #  
    MCC_v$DS_v[which(MCC_v$inputs=="ND")][i] <- round(ds$MCCsd[2], digits=2)
    MCC_v$DS_v[which(MCC_v$inputs=="D")][i] <- round(ds_diag$MCCsd[2], digits=2)
    #  
    MCC_v$PSG[which(MCC_v$inputs=="ND")][i] <- round(metapsg$MCC[2], digits=2)
    MCC_v$PSG[which(MCC_v$inputs=="D")][i] <- round(metapsg_diag$MCC[2], digits=2)
    #
    MCC_v$PSG_v[which(MCC_v$inputs=="ND")][i] <- round(metapsg$MCCsd[2], digits=2)
    MCC_v$PSG_v[which(MCC_v$inputs=="D")][i] <- round(metapsg_diag$MCCsd[2], digits=2)
    #    
    MCC_v$PSG_aw[which(MCC_v$inputs=="ND")][i] <- round(metapsg_adapt$MCC[2], digits=2)
    MCC_v$PSG_aw[which(MCC_v$inputs=="D")][i] <- round(metapsg_diag_adapt$MCC[2], digits=2)
    #
    MCC_v$PSG_aw_v[which(MCC_v$inputs=="ND")][i] <- round(metapsg_adapt$MCCsd[2], digits=2)
    MCC_v$PSG_aw_v[which(MCC_v$inputs=="D")][i] <- round(metapsg_diag_adapt$MCCsd[2], digits=2)
	#
	MCC_v$ASSET[which(MCC_v$inputs=="ND")][i] <- round(asset$MCC[1], digits=2)
	MCC_v$ASSET[which(MCC_v$inputs=="D")][i] <- " "#
	#
	MCC_v$ASSET_v[which(MCC_v$inputs=="ND")][i] <- round(asset$MCCsd[1], digits=2)
	MCC_v$ASSET_v[which(MCC_v$inputs=="D")][i] <- " "#	
	#	
	MCC_v$PLACO[which(MCC_v$inputs=="ND")][i] <- round(placo$MCC[1], digits=2)
	MCC_v$PLACO[which(MCC_v$inputs=="D")][i] <- " "#
	#
	MCC_v$PLACO_v[which(MCC_v$inputs=="ND")][i] <- round(placo$MCCsd[1], digits=2)
	MCC_v$PLACO_v[which(MCC_v$inputs=="D")][i] <- " "#		
	#
  }
  write.table(MCC_v, file=paste0(resdir,name_out,"_variable.txt"), quote = FALSE, sep = " &\ ")
}

# names_folders <- c("batch_data_diag_1/", "batch_data_diag_2/", "batch_data_diag_3/", "batch_data_diag_4/", "batch_data_diag_5/")
names_folders <- c("batch_data_cor25_1/", "batch_data_cor25_2/", "batch_data_cor25_3/", "batch_data_cor25_4/", "batch_data_cor25_5/")
# names_folders <- c("batch_data_cor50_1/", "batch_data_cor50_2/", "batch_data_cor50_3/", "batch_data_cor50_4/", "batch_data_cor50_5/")

#resdir="/home1/6_AMLAP/PE/METAPSG/SIMU/beta_from_CovX/results/simu_IGS_30/rep50/G5_SNR_03_beta_01_rep50/"
#resdir="/home1/6_AMLAP/PE/METAPSG/SIMU/beta_from_CovX/results/simu_IGS_30/rep50/G5_SNR_06_beta_015_rep50/"

resdir="/home1/6_AMLAP/PE/METAPSG/SIMU/beta_from_CovX/results/simu_IGS_30/rep50/G50_SNR_06_beta_015_rep50/"

# name_out <- "MCC_table_main_paper_cor25_beta15"
#name_out <- "MCC_table_main_paper_cor50_beta15"

name_out <- "MCC_table_paper_IGS30_cor25_beta15"

#
#
f_create_table_groupe(names_folders, name_out)
f_create_table_variable(names_folders, name_out)


