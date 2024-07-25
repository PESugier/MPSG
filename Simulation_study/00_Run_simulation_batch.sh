#===============================================================================
#           USAGE:  ./00_Run_simulation_batch.sh simulation_batch_1.ini
#
#     DESCRIPTION:  - global codes to run a simulation batch to test different methods adapted for pleiotropy: main purpose - to test Meta-PSG performances 
#                   - parameters need to be specified by the user in a separated .ini file
#
#===============================================================================
#          AUTHOR:  Pierre-Emmanuel Sugier
#         VERSION:  Pierre-Emmanuel Sugier
#         CREATED:  2024-03-15
#        REVISION:  2024-03-15
# WARNINGS N BUGS:  
#                   Is there some???


## REQUIRED INFORMATIONS
##########################################################################################################################################
start=$(date +%s.%N)

ref_file_parameters="$1"

source $ref_file_parameters

echo $nseed

cd ${scriptsimudir}"results"
mkdir ${name_batch}


## RUN BATCH FOR EACH METHOD
##########################################################################################################################################

# PLACO
################################################
Rscript ${scriptsimudir}/script_PLACO.R \
      --path_out ${outputdir} \
      --n_rep ${nseed} \
      --effect ${effect} \
      --SNR ${SNR} \
      --GS ${GS} \
      --IGS ${IGS} \
      --n_study ${K} \
      --n_indiv ${N} \
      --g ${G} \
      --mg ${mg} \
      --corr ${var_corr}
# ASSET
################################################
Rscript ${scriptsimudir}/script_ASSET.R \
      --path_out ${outputdir} \
      --n_rep ${nseed} \
      --effect ${effect} \
      --SNR ${SNR} \
      --GS ${GS} \
      --IGS ${IGS} \
      --n_study ${K} \
      --n_indiv ${N} \
      --g ${G} \
      --mg ${mg} \
      --corr ${var_corr}
# META-PSG
################################################
Rscript ${scriptsimudir}/script_MetaPSG_adapt.R \
      --path_metapsg ${metapsgdir} \
      --path_out ${outputdir} \
      --n_rep ${nseed} \
      --effect ${effect} \
      --SNR ${SNR} \
      --GS ${GS} \
      --IGS ${IGS} \
      --n_study ${K} \
      --n_indiv ${N} \
      --g ${G} \
      --mg ${mg} \
      --corr ${var_corr}
Rscript ${scriptsimudir}/script_MetaPSG_v2_adapt.R \
      --path_metapsg ${metapsgdir} \
      --path_out ${outputdir} \
      --n_rep ${nseed} \
      --effect ${effect} \
      --SNR ${SNR} \
      --GS ${GS} \
      --IGS ${IGS} \
      --n_study ${K} \
      --n_indiv ${N} \
      --g ${G} \
      --mg ${mg} \
      --corr ${var_corr}
# GCPBayes (DS)
################################################
Rscript ${scriptsimudir}/script_DS_short.R \
      --path_out ${outputdir} \
      --n_rep ${nseed} \
      --effect ${effect} \
      --SNR ${SNR} \
      --GS ${GS} \
      --IGS ${IGS} \
      --n_study ${K} \
      --n_indiv ${N} \
      --g ${G} \
      --mg ${mg} \
      --corr ${var_corr}
# META-PSG Diag
################################################  
Rscript ${scriptsimudir}/script_MetaPSG_diag_adapt.R \
      --path_metapsg ${metapsgdir} \
      --path_out ${outputdir} \
      --n_rep ${nseed} \
      --effect ${effect} \
      --SNR ${SNR} \
      --GS ${GS} \
      --IGS ${IGS} \
      --n_study ${K} \
      --n_indiv ${N} \
      --g ${G} \
      --mg ${mg} \
      --corr ${var_corr}
Rscript ${scriptsimudir}/script_MetaPSG_v2_diag_adapt.R \
      --path_metapsg ${metapsgdir} \
      --path_out ${outputdir} \
      --n_rep ${nseed} \
      --effect ${effect} \
      --SNR ${SNR} \
      --GS ${GS} \
      --IGS ${IGS} \
      --n_study ${K} \
      --n_indiv ${N} \
      --g ${G} \
      --mg ${mg} \
      --corr ${var_corr}
# GCPBayes (DS) Diag
################################################        
Rscript ${scriptsimudir}/script_DS_short_diag.R \
      --path_out ${outputdir} \
      --n_rep ${nseed} \
      --effect ${effect} \
      --SNR ${SNR} \
      --GS ${GS} \
      --IGS ${IGS} \
      --n_study ${K} \
      --n_indiv ${N} \
      --g ${G} \
      --mg ${mg} \
      --corr ${var_corr}
# MetaPSG without adaptive weights
Rscript ${scriptsimudir}/script_MetaPSG.R \
      --path_metapsg ${metapsgdir} \
      --path_out ${outputdir} \
      --n_rep ${nseed} \
      --effect ${effect} \
      --SNR ${SNR} \
      --GS ${GS} \
      --IGS ${IGS} \
      --n_study ${K} \
      --n_indiv ${N} \
      --g ${G} \
      --mg ${mg} \
      --corr ${var_corr}
Rscript ${scriptsimudir}/script_MetaPSG_diag.R \
      --path_metapsg ${metapsgdir} \
      --path_out ${outputdir} \
      --n_rep ${nseed} \
      --effect ${effect} \
      --SNR ${SNR} \
      --GS ${GS} \
      --IGS ${IGS} \
      --n_study ${K} \
      --n_indiv ${N} \
      --g ${G} \
      --mg ${mg} \
      --corr ${var_corr}
# END SCRIPT