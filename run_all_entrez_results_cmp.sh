#./run_all_entrez_results_cmp.sh

all_data=(
"GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal"
"GSE105318_ENCFF439QFU_DLD1"
#"GSE105465_ENCFF777DUA_Caki2_vs_GSE105235_ENCFF235TGH_G401 TCGAkich_norm_kich"
"GSE105566_ENCFF358MNA_Panc1"
"GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460"
"GSE106022_ENCFF614EKT_RPMI7951_vs_GSE105491_ENCFF458OWO_SKMEL5"
"GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP"
)

 script_name="all_entrez_fastSemSim_otherTADfile.R"
 
 for data in "${all_data[@]}"; do
    echo "> START for $data"
	echo Rscript $script_name $data
	Rscript $script_name $data
done

