*** MALARIA SNP ANALYSIS

*** SETTINGS AND PATHS
    clear all
	set more off
    capture log close
	program drop _all
	
	capture cd ""
	display _rc
	if _rc==0 {
		global jdlaptop = 1
		}
	else {
		global jdlaptop = 0
		}

	capture cd "" /* Run on MyResearch */
	display _rc

	local datapath  	= "Data/"
 	local logpath   	= "Log/"
	global outputpath 	= "Output/" /* Must be global to export in programs */
	
 ** Create date local to name files
	local currdate_raw: display %td_CCYY_NN_DD date(c(current_date), "DMY")
	local currdate = subinstr(trim("`currdate_raw'"), " " , "-", .)
		
	local inputdate = "`currdate'"
	
	local inputfile	 		= "Malaria_TAPBP_Variables_`inputdate'"
		
 ** Set exposure list

 ** Set local to change covariate list to include site_bantu instead of siteid for sensitivity analyses
	local ifsitebantu = 0

 ** File Name - General
	local filename = "TAPBP"
	
	global excelfile	  	= "Malaria_`filename'_Analysis_`currdate'.xlsx"  /* Must be global to export in programs */
	local logfile	  		= "Malaria_`filename'_Analysis_`currdate'"
	
 ** Standard covariate list
	local cov = "male eir_ln ib3.siteid" /* Different age var in each model */
	
 ** Site/Bantu combined variable covariate list
	if `ifsitebantu' == 1 {
		local cov = "male eir_ln ib3.site_bantu"/* Different age var in each model */
		global excelfile	  	= "Malaria_`filename'_Analysis_sitebantu_`currdate'.xlsx"  /* Must be global to export in programs */
		local logfile	  		= "Malaria_`filename'_Analysis_sitebantu_`currdate'"
	}
	
*** OPEN LOG FILE
	log using "`logpath'/`logfile'", replace
	
 ** Set max iterations
	set maxiter 500
	
********************************************************************************	
********************************************************************************

*** Call data
	use "`datapath'/`inputfile'", clear	

********************************************************************************	
	
*** Exposures
	* 3’: AA (high expression) vs. AG/GG (low expression)
	* 5’: CG/GG (high expression) vs. CC (low expression)
	* Binary dependence splits: 1.70 for ABC, 1.10 for A, 1.40 for B and 0.80 for C
		* Dependence = High dependence, Independence = Low dependence

 ** All exposures
	local exp_snp_exp 	exp_3_c snp_aa snp_cg_gg
	local exp_dep_c 	dep_all dep_hla_a dep_hla_b dep_hla_c 
	local exp_dep_3		i.dep_all_170_3_a i.dep_hla_a_110_3_a i.dep_hla_b_140_3_a i.dep_hla_c_080_3_a ///
						i.dep_all_170_3_c i.dep_hla_a_110_3_c i.dep_hla_b_140_3_c i.dep_hla_c_080_3_c
	local exp_dep_b 	dep_all_170_b dep_hla_a_110_b dep_hla_b_140_b dep_hla_c_080_b
	
	*local exp_snp_exp_m	exp_3_c_m snp_aa_m snp_cg_gg_m
	
	local exp_list_all = "snp_exp dep_c dep_b dep_3"
	
*** Subgroups
	local subgroups 	all_170 hla_a_110 hla_b_140 hla_c_080

********************************************************************************

*** Descriptive tables
	keep id hhid siteid gender child dep* exp* snp*
	
	duplicates drop
	
	tab snp_aa
	tab snp_cg_gg
	tab exp_3_c

	tab snp_aa siteid, col chi
	tab snp_cg_gg siteid, col chi
	tab exp_3_c siteid, col chi

	tab snp_aa child, col chi
	tab snp_cg_gg child, col chi
	tab exp_3_c child, col chi

	use "`datapath'/`inputfile'", clear	
	
********************************************************************************

*** Program to fill in exposure across the top of each sheet
	program label_reg
	if $jdlaptop == 1 {
		putexcel set "${outputpath}/${excelfile}", sheet("${outcome}_${explistname}") modify
		if $i == 1 {	
			putexcel C1 = "$exp"
			}
		else if $i > 1 {
			excelcol $colnum
			local colname `r(column)'
			putexcel `colname'1 = "$exp"
		}
	}
	end
 
*** Program to export regression results
	program reg_exp
		args b
		scalar econverged = e(converged) /* For reg_exp_one program below */
		
		mat mat1 = r(table)
		mat mat2 = mat1["b",1...]\mat1["pvalue",1...]\mat1["ll",1...]\mat1["ul",1...]
		mat mat3 = mat2'
		
	if $jdlaptop == 1 {
		putexcel set "${outputpath}/${excelfile}", sheet("${outcome}_${explistname}") modify	
		if $i == 1 {
			global colnum = 1 /* Start in first column */
			excelcol $colnum /* Convert column number to letter */
			local colname `r(column)' /* Create local of column letter */
			local rownum = $rownum + 1 /* Add 1 to row to get row 3 */
			putexcel `colname'$rownum = matrix(mat3), nformat("0.000") names /* Fill in matrix with column names */
			
			* Overwrite exposure name for that model with generic name for sheet
			putexcel B`rownum' = "Exposure"
			
			global colnum = $colnum + 2 /* Move two columns over from current column */
			excelcol $colnum  /* Convert column number to letter */
			local colname `r(column)' /* Create local of column letter */
			putexcel `colname'`=`rownum'-1' = "`b'" /* Fill in label for coef/HR/IRR/beta */
			local rownum = $rownum +`= rowsof(mat3)'+1 /* Add the length of the matrix + 1 to move down rows */
			putexcel B`rownum' = "N" /* Fill in labels for N, # clusters, # individuals */
			putexcel B`=`rownum'+1' = "# HH"
			putexcel B`=`rownum'+2' = "# ID"
		}
		
		else if $i > 1 {
			excelcol $colnum /* Start in current column */
			local colname `r(column)' /* Create local of column letter */
			putexcel `colname'$rownum = matrix(mat3), nformat("0.000") colnames /* Fill in matrix with ONLY column names (no row labels) */
			putexcel `colname'$rownum = "`b'" /* Fill in label for coef/HR/IRR/beta */
		}
			
		global rownum = $rownum +`= rowsof(mat3)'+1 /* Add the length of the matrix + 1 to move down rows */
		mat mat4 = [e(N)],[e(N_g)] /* Capture N and number of groups */
		mat mat5 = mat4' /* Transpose the matrix */
		putexcel `colname'$rownum = matrix(mat5) /* Fill in N and number of groups */
	}
	
	end
	
*** Program to export results in one row to merge with permutation results
	global j = 1 /* Set row counter */
	global site_id = 0 /* Set flag for site id */
	
	program reg_exp_one
	
		putexcel set "${outputpath}/${excelfile}", sheet("results_one") modify
		* Add row labels
		
		if $j == 1 {
			putexcel A1 = "outcome"
			putexcel B1 = "outcome_var"
			putexcel C1 = "exp"
			putexcel D1 = "coef"
			putexcel E1 = "p"
			putexcel F1 = "coef_ci_lb"
			putexcel G1 = "coef_ci_ub"
			putexcel H1 = "n_obs"
			putexcel I1 = "n_hh"
			putexcel J1 = "n_id"
			putexcel K1 = "among"
			putexcel L1 = "converged"
			putexcel M1 = "base"
			putexcel N1 = "aic"
			putexcel O1 = "bic"
			global j = $j + 1
			}
		
		di "$exp"
		
		if substr("$exp", 1, 2) == "i." {
			global exp_noi = regexr("$exp", "i\.", "")
			}
		else if substr("$exp", 1, 2) == "ib" {
			global exp_noi = regexr("$exp", "ib[0-9]\.", "")
			}
		else {
			global exp_noi = "$exp"
			}
			/* Remove i. from exposure to use string to match coefficient column names */
			
		di "$exp_noi"
		*mat mat1 = r(table) /* Use matrix from reg_exp; if don't run reg_exp first, need to unstar */
		
		local matnametest : colnames mat1 /* Save column names (coefficients) into a local */
		di "`matnametest'"
		di wordcount("`matnametest'")
		local matnametest: subinstr local matnametest "$exp_noi" "$exp_noi", all count(local c) 
			/* Use substring function to count number of coefficients with exposure name in them;
			   this equals the number of columns to extract from the matrix */
		di `c' /* Check number of columns to extract */
		
		if "$outcome_var" != "para_ln" {
			mat mat2 = mat1["b",1..`c']\mat1["pvalue",1..`c']\mat1["ll",1..`c']\mat1["ul",1..`c']
			mat mat3 = mat2'
		}
		
		if "$outcome_var" == "para_ln" {
			/* Exponeniate coefs for log parasite density for interpretability */
			mat mat2 = mat1["b",1..`c']\mat1["pvalue",1..`c']\mat1["ll",1..`c']\mat1["ul",1..`c']
			mat mat2a = mat2
			foreach i in 1 3 4 {
				/* Exponentiate rows for b, ll, ul */
					forvalues j = 1/`c' {
						mat mat2a[`i', `j'] = exp(mat2[`i', `j'])
					}
				}
			mat mat3 = mat2a'
		}
		
		putexcel A$j = "$outcome"
		putexcel B$j = matrix(mat3), rownames
			
		scalar n_obs = e(N)
		mat n = e(N_g)
		scalar n_hh = n[1,1]
		scalar n_id = n[1,2]
		
		putexcel H$j = n_obs
		putexcel I$j = n_hh
		putexcel J$j = n_id
		capture putexcel K$j = "$highlow"
		putexcel L$j = econverged
		
		estat ic
		mat mat6 = r(S)
		scalar aic = mat6[1,5]
		scalar bic = mat6[1,6]
		
		putexcel N$j = aic
		putexcel O$j = bic
		
			if `c' > 1 { /* Fill in 2nd row of coefficients if more than one coefficient */
				global j = $j + 1
				putexcel A$j = "$outcome"
				putexcel H$j = n_obs
				putexcel I$j = n_hh
				putexcel J$j = n_id
				capture putexcel K$j = "$highlow"
				putexcel L$j = econverged
				}
			
			if `c' > 2 { /* Fill in 3rd row of coefficients if more than one coefficient */
				global j = $j + 1
				putexcel A$j = "$outcome"
				putexcel H$j = n_obs
				putexcel I$j = n_hh
				putexcel J$j = n_id
				capture putexcel K$j = "$highlow"
				putexcel L$j = econverged
				}
				
			if `c' > 3 { /* Fill in 4th row of coefficients if more than one coefficient */
				global j = $j + 1
				putexcel A$j = "$outcome"
				putexcel H$j = n_obs
				putexcel I$j = n_hh
				putexcel J$j = n_id
				capture putexcel K$j = "$highlow"
				putexcel L$j = econverged
				
				}
		
			global j = $j + 1
	end
	
*** Program to run/export AIC/BIC
	program aic_bic
	
	if $jdlaptop==1 {
		putexcel set "${outputpath}/${excelfile}", sheet("${outcome}_${explistname}") modify
		estat ic
		local rownum = $rownum
		putexcel B`=`rownum'+3' = "AIC"
		putexcel B`=`rownum'+4' = "BIC"
		mat mat6 = r(S)
		scalar aic = mat6[1,5]
		scalar bic = mat6[1,6]
		
		excelcol $colnum
		local colname `r(column)'
		putexcel `colname'`=`rownum'+3' = aic, nformat("0.00")
		putexcel `colname'`=`rownum'+4' = bic, nformat("0.00")
		global rownum = `rownum'+6
	}
	end
		
********************************************************************************

*** Loop through exposures
	foreach exp_list in `exp_list_all' {
		
		global explistname = "`exp_list'"
	
********************************************************************************
********************************************************************************

*** Log parasite density: Routine visits
	global outcome = "para_ln_r"
	
 ** Call data
	use "`datapath'/`inputfile'", clear		

 ** Keep routine visits
	keep if visittype != 1
	
 ** Keep if parasite density > 0
	keep if parasitedensity>0 & parasitedensity<.
	
 ** Empty models
	putexcel set "${outputpath}/${excelfile}", sheet(Empty) modify
	reg para_ln
	estat ic
	mixed para_ln || id:
	estat ic
	mixed para_ln || hhid: || id:
	estat ic /* Gives you AIC and BIC */

 ** Linear models
	if $jdlaptop == 1 {
		putexcel set "${outputpath}/${excelfile}", sheet("${outcome}_${explistname}") modify
		putexcel A1 = "Log parasitemia: Routine visits"
	}
	
	global colnum = 1
	global i = 1
	
	foreach exp in `exp_`exp_list'' {
	
		global rownum = 2
		global exp = "`exp'"
		label_reg
		
		/* Crude mixed model */
		mixed para_ln `exp' || hhid: || id:, reml 
		reg_exp Coef
		aic_bic	
			
		/* Adjusted mixed model */
		mixed para_ln `exp' c.age##i.child `cov' || hhid: || id:, reml
		reg_exp Coef
		reg_exp_one
		aic_bic	
		
		
		global colnum = $colnum + 4
		global i = $i + 1	
		
		if "`exp'" == "snp_ag_gg" | "`exp'" == "snp_aa" |"`exp'" == "snp_cg_gg" |"`exp'" == "exp_t" |"`exp'" == "exp_3_c" | "`exp'" == "exp_6_c"  {
		
			foreach dep in `subgroups' {
				
				forvalues hl = 0/1 {
				
					* Fill in correct label for reg_exp (full model) spreadsheets
					global rownum = 2
					
					if `hl' == 0 {
						global highlow = "`dep'_l"
					}
					else if `hl' == 1 {
						global highlow = "`dep'_h"
					}
					
					putexcel set "${outputpath}/${excelfile}", sheet("${outcome}_${explistname}") modify
					excelcol $colnum
					local colname `r(column)'
					putexcel `r(column)'1 = "`exp'_${highlow}"
					
					mixed para_ln `exp' if dep_`dep'_b == `hl' || hhid: || id:, reml
					reg_exp Coef
					aic_bic	
												
					mixed para_ln `exp' c.age##i.child `cov' if dep_`dep'_b == `hl' || hhid: || id:, reml
					reg_exp Coef
					reg_exp_one
					aic_bic	
					
					
					global colnum = $colnum + 4
					global i = $i + 1	
					
				}
		
			}
			
			global highlow = ""
		}
		
	}

********************************************************************************

*** Log parasite density: Malaria care-seeking visits (non-routine)
	global outcome = "para_ln_m"
	
 ** Call data
	use "`datapath'/`inputfile'", clear		

 ** Keep non-routine visits
	keep if visittype == 1
	
 ** Keep if parasite density > 0
	keep if parasitedensity>0 & parasitedensity<.
	
 ** Drop if do NOT have malaria/do NOT have fever
	drop if febrile==0
	
 ** Empty models
	reg para_ln
	estat ic
	mixed para_ln || id:
	estat ic
	mixed para_ln || hhid: || id:
	estat ic
	
 ** Linear models
	if $jdlaptop == 1 {
		putexcel set "${outputpath}/${excelfile}", sheet("${outcome}_${explistname}") modify
		putexcel A1 = "Log parasitemia: Non-routine visits"
	}
	
	global colnum = 1
	global i = 1
	
	foreach exp in `exp_`exp_list'' {
		
		global rownum = 2
		global exp = "`exp'"
		label_reg
		
		mixed para_ln `exp' || hhid: || id:, reml
		reg_exp Coef
		aic_bic
		
		mixed para_ln `exp' c.age##i.child `cov' || hhid: || id:, reml
		reg_exp Coef
		reg_exp_one
		aic_bic
		
		
		global colnum = $colnum + 4
		global i = $i + 1	
		
		if "`exp'" == "snp_ag_gg" | "`exp'" == "snp_aa" |"`exp'" == "snp_cg_gg" |"`exp'" == "exp_t" |"`exp'" == "exp_3_c" | "`exp'" == "exp_6_c"  {
		
			foreach dep in `subgroups' {
				
				forvalues hl = 0/1 {
				
					* Fill in correct label for reg_exp (full model) spreadsheets
					global rownum = 2
					
					if `hl' == 0 {
						global highlow = "`dep'_l"
					}
					else if `hl' == 1 {
						global highlow = "`dep'_h"
					}
					
					putexcel set "${outputpath}/${excelfile}", sheet("${outcome}_${explistname}") modify
					excelcol $colnum
					local colname `r(column)'
					putexcel `r(column)'1 = "`exp'_${highlow}"
					
					mixed para_ln `exp' if dep_`dep'_b == `hl' || hhid: || id:, reml
					reg_exp Coef
					aic_bic	
												
					mixed para_ln `exp' c.age##i.child `cov' if dep_`dep'_b == `hl' || hhid: || id:, reml
					reg_exp Coef
					reg_exp_one
					aic_bic
					
					
					global colnum = $colnum + 4
					global i = $i + 1	
					
				}
		
			}
			
			global highlow = ""
		}
		
	}	

********************************************************************************		

*** Incident Malaria: Year
	global outcome = "inc_mal"

 ** Call data
	use "`datapath'/`inputfile'", clear	

 ** Keep annual variables
	keep id hhid site* *bantu* snp* exp* dep* year incidentmalaria_year pt_year age_y agecat_y male eir_ln g6pd_r hbs_r alphathal_r child
	duplicates drop /* One row per person per year */
	
	drop if pt_year==0 /* Drop if person-time = 0 */

 ** Create offset
	gen pt_days_log = log(pt_year/365.25) /* Log person-time */
 
 ** Empty models
	poisson incidentmalaria_year, offset(pt_days_log) irr
	estat ic
	mepoisson incidentmalaria_year, offset(pt_days_log) || id:, irr
	estat ic
	mepoisson incidentmalaria_year, offset(pt_days_log) || hhid: || id:, irr
	estat ic

 ** Poisson models
	if $jdlaptop == 1 {
		putexcel set "${outputpath}/${excelfile}", sheet("${outcome}_${explistname}") modify
		putexcel A1 = "Incident malaria: Annual"
	}
	
	global colnum = 1
	global i = 1
	
	foreach exp in `exp_`exp_list'' {
		
		global rownum = 2
		global exp = "`exp'"
		label_reg
		
		mepoisson incidentmalaria_year `exp', offset(pt_days_log) || hhid: || id:, irr
		reg_exp IRR
		aic_bic
		
		mepoisson incidentmalaria_year `exp' c.age_y##i.child `cov', offset(pt_days_log) || hhid: || id:, irr
		reg_exp IRR
		reg_exp_one
		aic_bic
		
		
		global colnum = $colnum + 4
		global i = $i + 1	
		
		if "`exp'" == "snp_ag_gg" | "`exp'" == "snp_aa" |"`exp'" == "snp_cg_gg" |"`exp'" == "exp_t" |"`exp'" == "exp_3_c" | "`exp'" == "exp_6_c"  {
		
			foreach dep in `subgroups' {
				
				forvalues hl = 0/1 {
				
					* Fill in correct label for reg_exp (full model) spreadsheets
					global rownum = 2
					
					if `hl' == 0 {
						global highlow = "`dep'_l"
					}
					else if `hl' == 1 {
						global highlow = "`dep'_h"
					}
					
					putexcel set "${outputpath}/${excelfile}", sheet("${outcome}_${explistname}") modify
					excelcol $colnum
					local colname `r(column)'
					putexcel `r(column)'1 = "`exp'_${highlow}"
					
					mepoisson incidentmalaria_year `exp' if dep_`dep'_b == `hl', offset(pt_days_log) || hhid: || id:, irr
					reg_exp IRR
					aic_bic	
					
					mepoisson incidentmalaria_year `exp' c.age_y##i.child `cov' if dep_`dep'_b == `hl', offset(pt_days_log) || hhid: || id:, irr
					reg_exp IRR
					reg_exp_one
					aic_bic	
					
					
					global colnum = $colnum + 4
					global i = $i + 1	
					
				}
		
			}
			
			global highlow = ""
		}
	
	}

********************************************************************************
********************************************************************************
	
*** Parasite Prevalence: Quarter
	global outcome = "para_prev"

 ** Call data
	use "`datapath'/`inputfile'", clear		
	
 ** Keep quarterly variables
	keep id hhid site* *bantu* snp* exp* dep* year q para_q age_q agecat_q male eir_ln g6pd_r hbs_r alphathal_r child
	duplicates drop /* One row per person per quarter */

 ** Empty models
	logistic para_q
	estat ic
	melogit para_q || id:, or
	estat ic
	melogit para_q || hhid: || id:, or
	estat ic
	
 ** Logistic models
	if $jdlaptop == 1 {
		putexcel set "${outputpath}/${excelfile}", sheet("${outcome}_${explistname}") modify
		putexcel A1 = "Parasite prevalence: Quarterly"
	}
	
	global colnum = 1
	global i = 1
	
	foreach exp in `exp_`exp_list'' {
		
		global rownum = 2
		global exp = "`exp'"
		label_reg
		
		melogit para_q `exp' || hhid: || id:, or
		reg_exp OR
		aic_bic
		
		melogit para_q `exp' c.age_q##i.child `cov' || hhid: || id:, or
		reg_exp OR
		reg_exp_one
		aic_bic
		
		
		global colnum = $colnum + 4
		global i = $i + 1	
		
		if "`exp'" == "snp_ag_gg" | "`exp'" == "snp_aa" |"`exp'" == "snp_cg_gg" |"`exp'" == "exp_t" |"`exp'" == "exp_3_c" | "`exp'" == "exp_6_c"  {
		
			foreach dep in `subgroups' {
				
				forvalues hl = 0/1 {
				
					* Fill in correct label for reg_exp (full model) spreadsheets
					global rownum = 2
					
					if `hl' == 0 {
						global highlow = "`dep'_l"
					}
					else if `hl' == 1 {
						global highlow = "`dep'_h"
					}
					
					putexcel set "${outputpath}/${excelfile}", sheet("${outcome}_${explistname}") modify
					excelcol $colnum
					local colname `r(column)'
					putexcel `r(column)'1 = "`exp'_${highlow}"
					
					melogit para_q `exp' if dep_`dep'_b == `hl' || hhid: || id:, or
					reg_exp OR
					aic_bic	
					
					melogit para_q `exp' c.age_q##i.child `cov' if dep_`dep'_b == `hl' || hhid: || id:, or
					reg_exp OR
					reg_exp_one
					aic_bic	
					
					
					global colnum = $colnum + 4
					global i = $i + 1	
					
				}
		
			}
			
			global highlow = ""
		}
	}
}
