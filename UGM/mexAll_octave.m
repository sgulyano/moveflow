mkoctfile --mex  -o minFunc_2009/lbfgsC.mex minFunc_2009/lbfgsC.c
mkoctfile --mex  -o minFunc_2009/mcholC.mex minFunc_2009/mcholC.c

mkoctfile --mex  -o KPM/max_mult.mex KPM/max_mult.c

mkoctfile --mex -Imex -o compiled/UGM_makeEdgeVEC.mex mex/UGM_makeEdgeVEC.c
mkoctfile --mex -Imex -o compiled/UGM_Decode_ExactC.mex mex/UGM_Decode_ExactC.c
mkoctfile --mex -Imex -o compiled/UGM_Infer_ExactC.mex mex/UGM_Infer_ExactC.c
mkoctfile --mex -Imex -o compiled/UGM_Infer_ChainC.mex mex/UGM_Infer_ChainC.c
mkoctfile --mex -Imex -o compiled/UGM_makeClampedPotentialsC.mex mex/UGM_makeClampedPotentialsC.c
mkoctfile --mex -Imex -o compiled/UGM_Decode_ICMC.mex mex/UGM_Decode_ICMC.c
mkoctfile --mex -Imex -o compiled/UGM_Decode_GraphCutC.mex mex/UGM_Decode_GraphCutC.c
mkoctfile --mex -Imex -o compiled/UGM_Sample_GibbsC.mex mex/UGM_Sample_GibbsC.c
mkoctfile --mex -Imex -o compiled/UGM_Infer_MFC.mex mex/UGM_Infer_MFC.c
mkoctfile --mex -Imex -o compiled/UGM_Infer_LBPC.mex mex/UGM_Infer_LBPC.c
mkoctfile --mex -Imex -o compiled/UGM_Decode_LBPC.mex mex/UGM_Decode_LBPC.c
mkoctfile --mex -Imex -o compiled/UGM_Infer_TRBPC.mex mex/UGM_Infer_TRBPC.c
mkoctfile --mex -Imex -o compiled/UGM_Decode_TRBPC.mex mex/UGM_Decode_TRBPC.c
mkoctfile --mex -Imex -o compiled/UGM_CRF_makePotentialsC.mex mex/UGM_CRF_makePotentialsC.c
mkoctfile --mex -Imex -o compiled/UGM_CRF_PseudoNLLC.mex mex/UGM_CRF_PseudoNLLC.c
mkoctfile --mex -Imex -o compiled/UGM_LogConfigurationPotentialC.mex mex/UGM_LogConfigurationPotentialC.c
mkoctfile --mex -Imex -o compiled/UGM_Decode_AlphaExpansionC.mex mex/UGM_Decode_AlphaExpansionC.c
mkoctfile --mex -Imex -o compiled/UGM_Decode_AlphaExpansionBetaShrinkC.mex mex/UGM_Decode_AlphaExpansionBetaShrinkC.c
mkoctfile --mex -Imex -o compiled/UGM_CRF_NLLC.mex mex/UGM_CRF_NLLC.c
mkoctfile --mex -Imex -o compiled/UGM_CRF_NLL_HiddenC.mex mex/UGM_CRF_NLL_HiddenC.c
