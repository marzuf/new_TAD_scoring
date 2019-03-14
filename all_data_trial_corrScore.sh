#!/usr/bin/bash

# ./all_data_trial_corrScore.sh


start_time=$(date -R)   


# Rscript allChromo_trials_corrScore.R 


scriptCoexpr="allChromo_trials_corrScore.R"

# should have hicds + expr ds !!!

all_TAD_files_ds=(
# DS launched 16.01.2018:
 "ENCSR079VIJ_G401_40kb TCGAkich_norm_kich"
 "ENCSR401TBQ_Caki2_40kb TCGAkich_norm_kich"
 "ENCSR079VIJ_G401ENCSR401TBQ_Caki2_40kb TCGAkich_norm_kich"

 "ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf"
 "ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutBRAF"
 "ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutCTNNB1"

 "ENCSR312KHQ_SK-MEL-5ENCSR862OGI_RPMI-7951_40kb TCGAskcm_lowInf_highInf"
 "ENCSR312KHQ_SK-MEL-5ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutBRAF"
 "ENCSR312KHQ_SK-MEL-5ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutCTNNB1"

 "ENCSR862OGI_RPMI-7951_40kb TCGAskcm_lowInf_highInf"
 "ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutBRAF"
 "ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutCTNNB1"

 "ENCSR346DCU_LNCaP_40kb TCGAprad_norm_prad"

 "ENCSR444WCZ_A549_40kb TCGAluad_mutKRAS_mutEGFR"
 "ENCSR444WCZ_A549_40kb TCGAluad_nonsmoker_smoker"
 "ENCSR444WCZ_A549_40kb TCGAluad_wt_mutKRAS"
 "ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc"

 "NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR"
 "NCI-H460_40kb TCGAluad_nonsmoker_smoker"
 "NCI-H460_40kb TCGAluad_wt_mutKRAS"
 "NCI-H460_40kb TCGAlusc_norm_lusc"

 "ENCSR444WCZ_A549NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR"
 "ENCSR444WCZ_A549NCI-H460_40kb TCGAluad_nonsmoker_smoker"
 "ENCSR444WCZ_A549NCI-H460_40kb TCGAluad_wt_mutKRAS"
 "ENCSR444WCZ_A549NCI-H460_40kb TCGAlusc_norm_lusc"

 "K562_40kb  TCGAlaml_wt_mutFLT3"

 "MCF-7_40kb TCGAbrca_lum_bas"

 "ENCSR549MGQ_T47D_40kb TCGAbrca_lum_bas"

 "MCF-7ENCSR549MGQ_T47D_40kb TCGAbrca_lum_bas"

 "Panc1_rep12_40kb TCGApaad_wt_mutKRAS"


# added 19.01.2018
 "GSE73782_PC3_40kb TCGAprad_norm_prad"
 "GSE73782_PC3_ICE_40kb TCGAprad_norm_prad"
 "GSE105318_DLD1_40kb TCGAcoad_msi_mss"
 "GSE105194_spinal_cord_40kb TCGAgbm_classical_mesenchymal"
 "GSE105194_spinal_cord_40kb TCGAgbm_classical_neural"
 "GSE105194_spinal_cord_40kb TCGAgbm_classical_proneural"
 "GSE105194_spinal_cord_40kb TCGAlgg_IDHwt_IDHmutnc"
 "ENCSR346DCU_LNCaPGSE73782_PC3_40kb TCGAprad_norm_prad"
 "ENCSR346DCU_LNCaPGSE73782_PC3_ICE_40kb TCGAprad_norm_prad"

# missing:
# colorectal; liver;  astrocytes 2/2; astrocyte consensus

# RUN 21.01.2018
 "GSE105194_cerebellum_40kb TCGAgbm_classical_mesenchymal"
 "GSE105194_cerebellum_40kb TCGAgbm_classical_neural"
 "GSE105194_cerebellum_40kb TCGAgbm_classical_proneural"
 "GSE105194_cerebellum_40kb TCGAlgg_IDHwt_IDHmutnc"
 
 "GSE105194_spinal_cordGSE105194_cerebellum_40kb TCGAgbm_classical_mesenchymal"
 "GSE105194_spinal_cordGSE105194_cerebellum_40kb TCGAgbm_classical_neural"
 "GSE105194_spinal_cordGSE105194_cerebellum_40kb TCGAgbm_classical_proneural"
 "GSE105194_spinal_cordGSE105194_cerebellum_40kb TCGAlgg_IDHwt_IDHmutnc"
 
 "GSE105381_HepG2_40kb TCGAlihc_wt_mutCTNNB1"

## added 01.03.2018 for variance comparison
"ENCSR444WCZ_A549_40kb TCGAluad_norm_luad"
"NCI-H460_40kb TCGAluad_norm_luad"
"ENCSR444WCZ_A549NCI-H460_40kb TCGAluad_norm_luad"
"GSE105381_HepG2_40kb TCGAlihc_norm_lihc"
)


# Rscript create_coexpr_sortNoDup_otherTADfile.R ENCSR346DCU_LNCaPGSE73782_PC3_ICE_40kb TCGAprad_norm_prad

for ds in "${all_TAD_files_ds[@]}"; do
    
    echo Rscript $scriptCoexpr $ds
    Rscript $scriptCoexpr $ds
done


###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################

echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0

