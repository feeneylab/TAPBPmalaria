*** MALARIA KIR/HLA ANALYSIS

*** SETTINGS AND PATHS
    clear
	clear matrix
	clear mata
	set more off
	set maxvar 20000
    capture log close
	program drop _all
	
	local path			= ""
	
	global datapath  	= "`path'/Data"
 	local logpath   	= "`path'/Log"
	global outputpath 	= "`path'/Output"
	
 ** Create date local to name files
	local currdate_raw: display %td_CCYY_NN_DD date(c(current_date), "DMY")
	local currdate = subinstr(trim("`currdate_raw'"), " " , "-", .)

	local inputfile_bp	 	= "tapbp_output_2021-03-04.csv"
	local inputfile_op	 	= "Malaria_TAPBP_Analysis_2021-03-02.xlsx"
	
	global outputfile_all	= "Malaria_TAPBP_Analysis_FDR_`currdate'"
	*local outputtables	 	= "Malaria_TAPBP_Tables_`currdate'"
	
	local logfile	  		= "Malaria_TAPBP_Analysis_FDR_`currdate'"

*** OPEN LOG FILE
	log using "`logpath'/`logfile'", replace

********************************************************************************

*** Create program to merge original and permuted models

program merge_o_p 
	args inputfile_op inputfile_bp outputfile outcomes sheetname

*** Call data: Original model output
	import excel "${outputpath}/`inputfile_op'", clear sheet("results_one") firstrow
	
 ** Drop reference group coefficients	
	gen base_drop = regexm(exp, "b\.") /* Tag base (reference) group coefficient */
	drop if base_drop == 1 /* Drop if reference group */
	drop base_drop

 ** Prepare for merge
  * Update with bases from exp list in analysis
	replace base=0 if regexm(exp, "dep.*_3_.") & base==. /* Fill in base = 0 for main models */
	
	tostring base, replace
	
forvalues i=0/2 {
		replace exp = exp + "_b" + base + "_`i'" if substr(exp, 1, 2)=="`i'."
		replace exp = substr(exp, 3, .) if substr(exp, 1, 2)=="`i'."
	}
	
	gen coef_name = outcome + "_" + exp if among==""
		replace coef_name = outcome + "_" + exp + "_" + among if among!=""	
	
	destring base, replace
	
	/* For model that is not converging, send results to missing */
	replace p = . if converged==0
	replace coef = . if converged==0
	replace coef_ci_lb = . if converged==0
	replace coef_ci_ub = . if converged==0

 ** Save data 
	save "${datapath}/`outputfile'", replace

********************************************************************************

*** Call data: Batch permutation output
	import delimited "${outputpath}/`inputfile_bp'", clear varnames(1)

 ** Rename
	rename * bp_*
	rename bp_outcome outcome
	rename bp_exp exp
	rename bp_coef coef_name

 ** Merge to original model output
	merge 1:1 coef_name using "${datapath}/`outputfile'"

	drop if _m == 1

 ** FDR of original p
  * Split p-values by outcome
  * Adjust using simes method (Simes [1986]; Benjamini and Hochberg [1995]; Benjamini and Yekutieli[2001, first method])
  	
	gen p_q = . /* Create empty variable for q-values to be combined into */
	
	foreach x in `outcomes' {
		gen p_`x' = p if outcome=="`x'"
		qqvalue p_`x', method(simes) qvalue(p_q_`x')
		count if p_q_`x'<0.05
		replace p_q = p_q_`x' if outcome=="`x'" /* Combine into one variable */
	}
	
	browse p bp_p p_* p_q
 	
	gen bp_q = . /* Create empty variable for q-values to be combined into */
	
	foreach x in `outcomes' {
		gen bp_p_`x' = bp_p if outcome=="`x'"
		replace bp_p_`x'=1/bp_repsnm if bp_p_`x'==0
		/* Compute q value assuming bp_p = 1/non-missing if bp_p==0, most conservative assumption */
		qqvalue bp_p_`x', method(simes) qvalue(bp_q_`x')
		count if bp_q_`x'<0.05	
		replace bp_q = bp_q_`x' if outcome=="`x'" /* Combine into one variable */
	}
	
	browse coef_name outcome p bp_p bp_q bp_*

 ** Label 3 level categorical variables
	gen coefficient = ""
		replace coefficient = "Dep + AG/GG" if regexm(coef_name, "a_b._1$")==1
		replace coefficient = "Dep + AA" if regexm(coef_name, "a_b._2$")==1
		replace coefficient = "Dep + CC" if regexm(coef_name, "c_b._1$")==1
		replace coefficient = "Dep + CG/GG" if regexm(coef_name, "c_b._2$")==1
	
	gen reference = "Independent" if base == 0 
		replace reference = "Dep + AG/GG" if regexm(coef_name, "a_b1_.$")==1
		replace reference = "Dep + CC" if regexm(coef_name, "c_b1_.$")==1
	
	gen note = "MODELS DO NOT CONVERGE" if coef == .
	
 ** Save data 
	sort outcome exp

	export excel coef_name outcome exp among coefficient reference coef coef_ci_lb coef_ci_ub p bp_p bp_q n_id note using "${outputpath}/${outputfile_all}.xlsx", sheetreplace sheet(`sheetname') firstrow(variables)
	
end 

********************************************************************************

*** Main outcomes

*** Merge data
	merge_o_p ///
	/*inputfile_op*/	"`inputfile_op'" ///
	/*inputfile_bp*/	"`inputfile_bp'" ///
	/*outputfile*/ 		"Malaria_TAPBP_Analysis_FDR_`currdate'" ///
	/*outcomes*/ 		"inc_mal para_ln_m para_ln_r para_prev" ///
	/*sheetname*/ 		"main"
