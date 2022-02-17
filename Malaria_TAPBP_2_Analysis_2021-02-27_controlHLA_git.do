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
	
 ** File Name - Primary
	local filename = "TAPBP"
	
	global excelfile	  	= "Malaria_`filename'_Analysis_`currdate'_controlHLA.xlsx"  /* Must be global to export in programs */
	local logfile	  		= "Malaria_`filename'_Analysis_`currdate'_controlHLA"
	
 ** Standard covariate list
	local cov = "male eir_ln ib3.siteid" /* Different age var in each model */
	
 ** Site/Bantu combined variable covariate list
  * Set local to change covariate list to include site_bantu instead of siteid for sensitivity analyses
	local ifsitebantu = 0
	
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
	
*** Program to export results in one row to merge with permutation results
	global j = 1 /* Set row counter */
	global site_id = 0 /* Set flag for site id */
	
	program reg_exp_one
	
		scalar econverged = e(converged) /* For reg_exp_one program below */
	
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
			putexcel P1 = "B5301"
			putexcel Q1 = "C0602"
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
		mat mat1 = r(table) /* Use matrix from reg_exp; if don't run reg_exp first, need to unstar */
		
		local matnametest : colnames mat1 /* Save column names (coefficients) into a local */
		di "`matnametest'"
		di wordcount("`matnametest'")
		local matnametest: subinstr local matnametest "$exp_noi" "$exp_noi", all count(local c) 
			/* Use subinstring function to count number of coefficients with exposure name in them;
			   this equals the number of columns to extract from the matrix */
		di `c' /* Check number of columns to extract */
		
		local matnametest: subinstr local matnametest "hla_b_l_5301" "hla_b_l_5301", all count(local c_b5301) 
			/* Use subinstring function to identify whether B5301 was included in model */
		di `c_b5301' /* Check whether it identifies B5301 in covariate list */
		
		local matnametest: subinstr local matnametest "hla_c_l_0602" "hla_c_l_0602", all count(local c_c0602) 
			/* Use subinstring function to identify whether C0602 was included in model */
		di `c_c0602' /* Check whether it identifies C0602 in covariate list */
		
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
		
		if `c_b5301'==1 {
			putexcel P$j = 1
		}
		
		if `c_c0602'==1 {
			putexcel Q$j = 1
		}
		
			if `c' > 1 { /* Fill in 2nd row of coefficients if more than one coefficient */
				global j = $j + 1
				putexcel A$j = "$outcome"
				putexcel H$j = n_obs
				putexcel I$j = n_hh
				putexcel J$j = n_id
				capture putexcel K$j = "$highlow"
				putexcel L$j = econverged
				putexcel N$j = aic
				putexcel O$j = bic
				if `c_b5301'==1 {
					putexcel P$j = 1
				}
				
				if `c_c0602'==1 {
					putexcel Q$j = 1
				}
			}
			
			if `c' > 2 { /* Fill in 3rd row of coefficients if more than one coefficient */
				global j = $j + 1
				putexcel A$j = "$outcome"
				putexcel H$j = n_obs
				putexcel I$j = n_hh
				putexcel J$j = n_id
				capture putexcel K$j = "$highlow"
						putexcel L$j = econverged
				putexcel N$j = aic
				putexcel O$j = bic
				if `c_b5301'==1 {
					putexcel P$j = 1
				}
				
				if `c_c0602'==1 {
					putexcel Q$j = 1
				}
			}
				
			if `c' > 3 { /* Fill in 4th row of coefficients if more than one coefficient */
				global j = $j + 1
				putexcel A$j = "$outcome"
				putexcel H$j = n_obs
				putexcel I$j = n_hh
				putexcel J$j = n_id
				capture putexcel K$j = "$highlow"
						putexcel L$j = econverged
				putexcel N$j = aic
				putexcel O$j = bic
				if `c_b5301'==1 {
					putexcel P$j = 1
				}
				
				if `c_c0602'==1 {
					putexcel Q$j = 1
				}
			}
		
			global j = $j + 1
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
	reg para_ln
	estat ic
	mixed para_ln || id:
	estat ic
	mixed para_ln || hhid: || id:
	estat ic /* Gives you AIC and BIC */

 ** Linear models
	foreach exp in `exp_`exp_list'' {

		global exp = "`exp'"
		
		/* Crude mixed model - all variables */
		mixed para_ln `exp' || hhid: || id:, reml 	
			
		/* Adjusted models */
		if strmatch("`exp'", "*_all*")==0 & strmatch("`exp'", "*hla_b*")==0 & strmatch("`exp'", "*hla_c*")==0 {
			/* Adjusted mixed model  - NOT CONTROLLING FOR HLA */
			mixed para_ln `exp' c.age##i.child `cov' || hhid: || id:, reml
			reg_exp_one
		}
		
		if strmatch("`exp'", "*_all*")==1 {
			/* Adjusted mixed model  - CONTROLLING FOR BOTH HLA */
			mixed para_ln `exp' c.age##i.child `cov' hla_b_l_5301 hla_c_l_0602 || hhid: || id:, reml
			reg_exp_one	
		}
		
		if strmatch("`exp'", "*hla_b*")==1 {
			/* Adjusted mixed model  - CONTROLLING FOR HLA-B*53:01 */
			mixed para_ln `exp' c.age##i.child `cov' hla_b_l_5301 || hhid: || id:, reml
			reg_exp_one	
		}
		
		if strmatch("`exp'", "*hla_c*")==1 {
			/* Adjusted mixed model  - CONTROLLING FOR HLA-C*06:02 */
			mixed para_ln `exp' c.age##i.child `cov' hla_c_l_0602 || hhid: || id:, reml
			reg_exp_one	
		}
			
		
		/* Subgroup models */ 
		if "`exp'" == "snp_aa"|"`exp'" == "snp_cg_gg" |"`exp'" == "exp_3_c" {
			foreach dep in `subgroups' {
				
				di "`dep'"
				
				forvalues hl = 0/1 {
				
					if `hl' == 0 {
						global highlow = "`dep'_l"
						di "$highlow"
					}
					else if `hl' == 1 {
						global highlow = "`dep'_h"
						di "$highlow"
					}
					
					mixed para_ln `exp' if dep_`dep'_b == `hl' || hhid: || id:, reml
										
					if strmatch("`dep'", "all*")==0 & strmatch("`dep'", "hla_b*")==0 & strmatch("`dep'", "hla_c*")==0 {
						/* Adjusted mixed model  - NOT CONTROLLING FOR HLA */
						mixed para_ln `exp' c.age##i.child `cov' if dep_`dep'_b == `hl' || hhid: || id:, reml
						reg_exp_one	
					}
					
					if strmatch("`dep'", "all*")==1 {
						/* Adjusted mixed model  - CONTROLLING FOR BOTH HLA */
						mixed para_ln `exp' c.age##i.child `cov' hla_b_l_5301 hla_c_l_0602 if dep_`dep'_b == `hl' || hhid: || id:, reml
						reg_exp_one		
					}
					
					if strmatch("`dep'", "hla_b*")==1 {
						/* Adjusted mixed model  - CONTROLLING FOR HLA-B*53:01 */
						mixed para_ln `exp' c.age##i.child `cov' hla_b_l_5301 if dep_`dep'_b == `hl' || hhid: || id:, reml
						reg_exp_one		
					}
					
					if strmatch("`dep'", "hla_c*")==1 {
						/* Adjusted mixed model  - CONTROLLING FOR HLA-C*06:02 */
						mixed para_ln `exp' c.age##i.child `cov' hla_c_l_0602 if dep_`dep'_b == `hl' || hhid: || id:, reml
						reg_exp_one	
					}
									
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
	estat ic /* Gives you AIC and BIC */

 ** Linear models
	foreach exp in `exp_`exp_list'' {

		global exp = "`exp'"
		
		/* Crude mixed model - all variables */
		mixed para_ln `exp' || hhid: || id:, reml 	
			
		/* Adjusted models */
		if strmatch("`exp'", "*_all*")==0 & strmatch("`exp'", "*hla_b*")==0 & strmatch("`exp'", "*hla_c*")==0 {
			/* Adjusted mixed model  - NOT CONTROLLING FOR HLA */
			mixed para_ln `exp' c.age##i.child `cov' || hhid: || id:, reml
			reg_exp_one
		}
		
		if strmatch("`exp'", "*_all*")==1 {
			/* Adjusted mixed model  - CONTROLLING FOR BOTH HLA */
			mixed para_ln `exp' c.age##i.child `cov' hla_b_l_5301 hla_c_l_0602 || hhid: || id:, reml
			reg_exp_one	
		}
		
		if strmatch("`exp'", "*hla_b*")==1 {
			/* Adjusted mixed model  - CONTROLLING FOR HLA-B*53:01 */
			mixed para_ln `exp' c.age##i.child `cov' hla_b_l_5301 || hhid: || id:, reml
			reg_exp_one	
		}
		
		if strmatch("`exp'", "*hla_c*")==1 {
			/* Adjusted mixed model  - CONTROLLING FOR HLA-C*06:02 */
			mixed para_ln `exp' c.age##i.child `cov' hla_c_l_0602 || hhid: || id:, reml
			reg_exp_one	
		}
			
		
		/* Subgroup models */ 
		if "`exp'" == "snp_aa"|"`exp'" == "snp_cg_gg" |"`exp'" == "exp_3_c" {
			foreach dep in `subgroups' {
				
				di "`dep'"
				
				forvalues hl = 0/1 {
				
					if `hl' == 0 {
						global highlow = "`dep'_l"
						di "$highlow"
					}
					else if `hl' == 1 {
						global highlow = "`dep'_h"
						di "$highlow"
					}
					
					mixed para_ln `exp' if dep_`dep'_b == `hl' || hhid: || id:, reml
										
					if strmatch("`dep'", "all*")==0 & strmatch("`dep'", "hla_b*")==0 & strmatch("`dep'", "hla_c*")==0 {
						/* Adjusted mixed model  - NOT CONTROLLING FOR HLA */
						mixed para_ln `exp' c.age##i.child `cov' if dep_`dep'_b == `hl' || hhid: || id:, reml
						reg_exp_one	
					}
					
					if strmatch("`dep'", "all*")==1 {
						/* Adjusted mixed model  - CONTROLLING FOR BOTH HLA */
						mixed para_ln `exp' c.age##i.child `cov' hla_b_l_5301 hla_c_l_0602 if dep_`dep'_b == `hl' || hhid: || id:, reml
						reg_exp_one		
					}
					
					if strmatch("`dep'", "hla_b*")==1 {
						/* Adjusted mixed model  - CONTROLLING FOR HLA-B*53:01 */
						mixed para_ln `exp' c.age##i.child `cov' hla_b_l_5301 if dep_`dep'_b == `hl' || hhid: || id:, reml
						reg_exp_one		
					}
					
					if strmatch("`dep'", "hla_c*")==1 {
						/* Adjusted mixed model  - CONTROLLING FOR HLA-C*06:02 */
						mixed para_ln `exp' c.age##i.child `cov' hla_c_l_0602 if dep_`dep'_b == `hl' || hhid: || id:, reml
						reg_exp_one	
					}
									
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
	keep id hhid site* *bantu* snp* exp* dep* year incidentmalaria_year pt_year age_y agecat_y male eir_ln g6pd_r hbs_r alphathal_r child hla*
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
	foreach exp in `exp_`exp_list'' {
		global exp = "`exp'"
		
		/* Crude mixed model - all variables */
		mepoisson incidentmalaria_year `exp', offset(pt_days_log) || hhid: || id:, irr
		
		
		/* Adjusted models */
		if strmatch("`exp'", "*_all*")==0 & strmatch("`exp'", "*hla_b*")==0 & strmatch("`exp'", "*hla_c*")==0 {
			/* Adjusted mixed model  - NOT CONTROLLING FOR HLA */
			mepoisson incidentmalaria_year `exp' c.age_y##i.child `cov', offset(pt_days_log) || hhid: || id:, irr
			reg_exp_one
		}
		
		if strmatch("`exp'", "*_all*")==1 {
			/* Adjusted mixed model  - CONTROLLING FOR BOTH HLA */
			mepoisson incidentmalaria_year `exp' c.age_y##i.child `cov' hla_b_l_5301 hla_c_l_0602, offset(pt_days_log) || hhid: || id:, irr
			reg_exp_one	
		}
		
		if strmatch("`exp'", "*hla_b*")==1 {
			/* Adjusted mixed model  - CONTROLLING FOR HLA-B*53:01 */
			mepoisson incidentmalaria_year `exp' c.age_y##i.child `cov' hla_b_l_5301, offset(pt_days_log) || hhid: || id:, irr
			reg_exp_one	
		}
		
		if strmatch("`exp'", "*hla_c*")==1 {
			/* Adjusted mixed model  - CONTROLLING FOR HLA-C*06:02 */
			mepoisson incidentmalaria_year `exp' c.age_y##i.child `cov' hla_c_l_0602, offset(pt_days_log) || hhid: || id:, irr
			reg_exp_one	
		}
			
		
		/* Subgroup models */ 
		if "`exp'" == "snp_aa"|"`exp'" == "snp_cg_gg" |"`exp'" == "exp_3_c" {
			foreach dep in `subgroups' {
				
				di "`dep'"
				
				forvalues hl = 0/1 {

					if `hl' == 0 {
						global highlow = "`dep'_l"
						di "$highlow"
					}
					else if `hl' == 1 {
						global highlow = "`dep'_h"
						di "$highlow"
					}
										
					mepoisson incidentmalaria_year `exp' if dep_`dep'_b == `hl', offset(pt_days_log) || hhid: || id:, irr
					
						
					/* Adjusted models */
					if strmatch("`dep'", "all*")==0 & strmatch("`dep'", "hla_b*")==0 & strmatch("`dep'", "hla_c*")==0 {
						/* Adjusted mixed model  - NOT CONTROLLING FOR HLA */
						mepoisson incidentmalaria_year `exp' c.age_y##i.child `cov' if dep_`dep'_b == `hl', offset(pt_days_log) || hhid: || id:, irr
						reg_exp_one
					}
					
					if strmatch("`dep'", "all*")==1 {
						/* Adjusted mixed model  - CONTROLLING FOR BOTH HLA */
						mepoisson incidentmalaria_year `exp' c.age_y##i.child `cov' hla_b_l_5301 hla_c_l_0602 if dep_`dep'_b == `hl', offset(pt_days_log) || hhid: || id:, irr
						reg_exp_one	
					}
					
					if strmatch("`dep'", "hla_b*")==1 {
						/* Adjusted mixed model  - CONTROLLING FOR HLA-B*53:01 */
						mepoisson incidentmalaria_year `exp' c.age_y##i.child `cov' hla_b_l_5301 if dep_`dep'_b == `hl', offset(pt_days_log) || hhid: || id:, irr
						reg_exp_one	
					}
					
					if strmatch("`dep'", "hla_c*")==1 {
						/* Adjusted mixed model  - CONTROLLING FOR HLA-C*06:02 */
						mepoisson incidentmalaria_year `exp' c.age_y##i.child `cov' hla_c_l_0602 if dep_`dep'_b == `hl', offset(pt_days_log) || hhid: || id:, irr
						reg_exp_one	
					}		
					
				}
		
			}
			
			global highlow = ""
		}
	
	}

********************************************************************************

*** Parasite Prevalence: Quarter
	global outcome = "para_prev"

 ** Call data
	use "`datapath'/`inputfile'", clear		
	
 ** Keep quarterly variables
	keep id hhid site* *bantu* snp* exp* dep* year q para_q age_q agecat_q male eir_ln g6pd_r hbs_r alphathal_r child hla*
	duplicates drop /* One row per person per quarter */
	
 ** Empty models
	logistic para_q
	estat ic
	melogit para_q || id:, or
	estat ic
	melogit para_q || hhid: || id:, or
	estat ic
	
 ** Logistic models
	foreach exp in `exp_`exp_list'' {
		
		global exp = "`exp'"

		/* Crude mixed model - all variables */
		melogit para_q `exp' || hhid: || id:, or
		
		/* Adjusted models */
		if strmatch("`exp'", "*_all*")==0 & strmatch("`exp'", "*hla_b*")==0 & strmatch("`exp'", "*hla_c*")==0 {
			/* Adjusted mixed model  - NOT CONTROLLING FOR HLA */
			melogit para_q `exp' c.age_q##i.child `cov' || hhid: || id:, or
			reg_exp_one
		}
		
		if strmatch("`exp'", "*_all*")==1 {
			/* Adjusted mixed model  - CONTROLLING FOR BOTH HLA */
			melogit para_q `exp' c.age_q##i.child `cov' hla_b_l_5301 hla_c_l_0602 || hhid: || id:, or
			reg_exp_one	
		}
		
		if strmatch("`exp'", "*hla_b*")==1 {
			/* Adjusted mixed model  - CONTROLLING FOR HLA-B*53:01 */
			melogit para_q `exp' c.age_q##i.child `cov' hla_b_l_5301 || hhid: || id:, or
			reg_exp_one	
		}
		
		if strmatch("`exp'", "*hla_c*")==1 {
			/* Adjusted mixed model  - CONTROLLING FOR HLA-C*06:02 */
			melogit para_q `exp' c.age_q##i.child `cov' hla_c_l_0602 || hhid: || id:, or
			reg_exp_one	
		}
			
		
		/* Subgroup models */ 
		if "`exp'" == "snp_aa"|"`exp'" == "snp_cg_gg" |"`exp'" == "exp_3_c" {
			
			foreach dep in `subgroups' {
				
				di "`dep'"
				
				forvalues hl = 0/1 {

					if `hl' == 0 {
						global highlow = "`dep'_l"
						di "$highlow"
					}
					else if `hl' == 1 {
						global highlow = "`dep'_h"
						di "$highlow"
					}
					
					
					melogit para_q `exp' if dep_`dep'_b == `hl' || hhid: || id:, or
						
					/* Adjusted models */
					if strmatch("`dep'", "all*")==0 & strmatch("`dep'", "hla_b*")==0 & strmatch("`dep'", "hla_c*")==0 {
						/* Adjusted mixed model  - NOT CONTROLLING FOR HLA */
						melogit para_q `exp' c.age_q##i.child `cov' if dep_`dep'_b == `hl' || hhid: || id:, or
						reg_exp_one
					}
					
					if strmatch("`dep'", "all*")==1 {
						/* Adjusted mixed model  - CONTROLLING FOR BOTH HLA */
						melogit para_q `exp' c.age_q##i.child `cov' hla_b_l_5301 hla_c_l_0602 if dep_`dep'_b == `hl' || hhid: || id:, or
						reg_exp_one	
					}
					
					if strmatch("`dep'", "hla_b*")==1 {
						/* Adjusted mixed model  - CONTROLLING FOR HLA-B*53:01 */
						melogit para_q `exp' c.age_q##i.child `cov' hla_b_l_5301 if dep_`dep'_b == `hl' || hhid: || id:, or
						reg_exp_one	
					}
					
					if strmatch("`dep'", "hla_c*")==1 {
						/* Adjusted mixed model  - CONTROLLING FOR HLA-C*06:02 */
						melogit para_q `exp' c.age_q##i.child `cov' hla_c_l_0602 if dep_`dep'_b == `hl' || hhid: || id:, or
						reg_exp_one	
					}
						
				}
		
			}
			
			global highlow = ""
		}
	}
} /* Close loop through exposure lists */
