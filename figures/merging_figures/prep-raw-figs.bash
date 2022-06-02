#/bin/bash 
resScale=150

# # folder abbreviations
# hvf=SW_HV_only
# raf=SW_Ray_phV_only 
# lof=SW_Lov_phV_only
# spf=RF_Sp_ccp_only 
# hkf=HKstack_P_only

# data abreviations
figmod=final_model 
figdat=final_true_vs_pred_data_wavs

#!/bin/bash
for fold in SW_HV_only SW_Ray_phV_only SW_Lov_phV_only RF_Sp_ccp_only HKstack_P_only
do
    mkdir -p ./merged_figs/$fold
    magick -density $(( $resScale * 1 )) $fold/final_model.pdf \
        ./merged_figs/$fold/final_model.jpg
    magick -density $(( $resScale * 1 )) $fold/final_true_vs_pred_data_wavs.pdf \
        ./merged_figs/$fold/final_true_vs_pred_data_wavs.jpg
done


# magick -density $(( $resScale * 1 )) SW_HV_only/ figs-prepped/final_model.jpg
# magick -density $(( $resScale * 1 )) figs-raw/prior2posterior.pdf figs-prepped/prior2posterior.jpg
