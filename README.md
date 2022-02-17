# TAPBPmalaria
Imputed TAPBP expression levels impact the course of malaria in an HLA allotype-dependent manner

Stata do files:

Malaria_TAPBP_1_Variables_2021-02-27_git.do Creates all variables necessary for analysis. Creates below .dta file.

Malaria_TAPBP_2_Analysis_2022-01-25_primaryanalysis_git.do Primary analysis file.

Malaria_TAPBP_2_Analysis_2021-02-27_controlHLA_git.do Same as above primary analysis but also controls for HLA-B*53:01 and HLA-C*06:02 for relevant models.

Malaria_TAPBP_3_Analysis_FDR_2021-03-08_git.do Merges primary analysis and batch permutation results. Calculates false discovery rate q-values.

BatchPermute_allmodels_TAPBP_2020-06-24_git.do Permutation analysis file.

Stata dta files:

Malaria_TAPBP_Variables_2022-02-16_git.dta Study data.
