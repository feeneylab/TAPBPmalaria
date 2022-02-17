*** Location of ritest package
	sysdir set PLUS ""

*** Set locals
	local exp : env varExp
	local envseed : env varSeed
	local outcome : env varOut
	local base : env varBase
	
	local hl : env varHL
	local dep : env varDep
	
	if "`hl'" == "0" {
		local hl_name "l"
		}
	if "`hl'" == "1" {
		local hl_name "h"
		}
	
*** Log/Output Settings
	capture log close
	
	if "`hl'"=="" & "`dep'"=="" {
		
		if "`base'"=="" { // Note base only applies to 3 level categorical models, never those subset by a binary dependence variable
			local filename = "`outcome'_`exp'"
		}
	
		else if "`base'"!="" {
			local filename = "`outcome'_`exp'_b`base'"
		}
	}
	
	else if "`hl'"!="" & "`dep'"!="" {
		local filename = "`outcome'_`exp'_`dep'_`hl_name'"
	}
	
	log using "Log/`filename'", replace
	local outcome_data	= "Output_Data/`filename'"
	local outcome_csv	= "Output/`filename'"
	
*** Settings
	di "`exp'"
	di "`base'"
	di "`outcome'"
	di "`dep'"
	di "`hl'"
	di "`hl_name'"
	di `envseed'
	
	local reps = 10000
	
	local seed = `envseed' + 80051
	set seed `seed'
	di `seed'
	set maxiter 500
	
*** Call Data
	local data_date = ""
	use "Data/Malaria_TAPBP_Variables_`data_date'.dta", clear		
	
*** Keep subset of data
	if "`hl'"!="" & "`dep'"!="" {
		keep if dep_`dep'_b == `hl'
		tab dep_`dep'_b
	}	
	
*** Standard covariate list
	local cov = "male eir_ln ib3.siteid" /* Different age var in each model */	
	
*** Set model variable and coefficient locals
	quietly tab `exp'
	
  * Binary/Continuous
	local exp_m = "`exp'"
	
	if "`hl'"=="" & "`dep'"=="" & regexm("`exp'", "dep.*_3_.") == 0 {
		local extract = "`outcome'_`exp' = _b[`exp']"
		di "`extract'"
		
	}
		
	if "`hl'"!="" & "`dep'"!="" & regexm("`exp'", "dep.*_3_.") == 0 {
		local extract = "`outcome'_`exp'_`dep'_`hl_name' = _b[`exp']"
		di "`extract'"
	}
	
  * Define categorical variables

	if regexm("`exp'", "dep.*_3_.") == 1 { /* For 3 level categorical dep/exp variables */
		
		levelsof `exp', local(levels)
			
		if "`base'" == "" {
			local exp_m = "i.`exp'" 
			di "`exp_m'"
			quietly sum `exp'
			local base = r(min)
			di `base'
			}	
		
		else {
			sum `exp'
			local exp_m = "ib`base'.`exp'"
			di "`exp_m'"
			}	
			
		local levels_extract: list levels- base /* List all levels except base */
		di "`levels_extract'"

		foreach l of local levels_extract { /* Loop through all levels except base to specify what to export */
			local extract_raw "`outcome'_`exp'_b`base'_`l' = _b[`l'.`exp']"
			
			di "`extract_raw'"
			local extract "`extract'" "`extract_raw' "
			di "`extract'"
			}
		}
		
********************************************************************************	
	
*** Log parasite density: Malaria visits

	if "`outcome'" == "para_ln_m" {

	*** Prepare data
	 ** Keep non-routine visits
		keep if visittype == 1
		
	 ** Keep if parasite density > 0
		keep if parasitedensity>0 & parasitedensity<.
		
	 ** Drop if do NOT have malaria/do NOT have fever
		drop if febrile==0

*** Permutations/Model		
	ritest `exp' "`extract'", reps(`reps') cluster(id) noisily reject(e(converged)!=1) saving("`outcome_data'", replace): ///
	mixed para_ln `exp_m' c.age##i.child `cov' || hhid: || id:, reml
		
	}
		
********************************************************************************	
		
*** Log parasite density: Routine visits	
	
	else if "`outcome'" == "para_ln_r" {
	
	*** Prepare data
	 ** Keep routine visits
		keep if visittype != 1
		
	 ** Keep if parasite density > 0
		keep if parasitedensity>0 & parasitedensity<.

*** Permutations/Model
	ritest `exp' "`extract'", reps(`reps') cluster(id) noisily reject(e(converged)!=1) saving("`outcome_data'", replace): ///
	mixed para_ln `exp_m' c.age##i.child `cov' || hhid: || id:, reml		
	
	}
	
********************************************************************************		
	
*** Incident Malaria: Yearly	

	else if "`outcome'" == "inc_mal" {

	*** Prepare data
	 ** Keep annual variables
		keep id hhid site* *bantu* snp* exp* dep* year incidentmalaria_year pt_year age_y agecat_y male eir_ln g6pd_r hbs_r alphathal_r child
		duplicates drop
		
		drop if pt_year==0

	 ** Create offset
		gen pt_days_log = log(pt_year/365.25)

*** Permutations/Model
	ritest `exp' "`extract'", reps(`reps') cluster(id) noisily reject(e(converged)!=1) saving("`outcome_data'", replace): ///
	mepoisson incidentmalaria_year `exp_m' c.age_y##i.child `cov', offset(pt_days_log) || hhid: || id:, irr
		
	}
		
********************************************************************************	

*** Parasite Prevalence: Quarter	
	
	else if "`outcome'" == "para_prev" {
	
	*** Prepare data
	 ** Keep quarterly variables
		keep id hhid site* *bantu* snp* exp* dep* year q para_q age_q agecat_q male eir_ln g6pd_r hbs_r alphathal_r child
		duplicates drop

*** Permutations/Model
	ritest `exp' "`extract'", reps(`reps') cluster(id) noisily reject(e(converged)!=1) saving("`outcome_data'", replace): ///
	melogit para_q `exp_m' c.age_q##i.child `cov' || hhid: || id:, or		

	}
	
********************************************************************************
	
*** Fill in data for export

	gen outcome = ""
	gen exp = ""
	gen coef = ""
	gen nobs = .
	gen obs = .
	gen cttrue = .
	gen repsn = .
	gen repsnm = .
	gen p = .
	gen cilb = .
	gen ciub = .
	gen seed = .
	gen high_low = ""
	gen dep = ""
	gen base = .
	
	local num_extract = r(k_exp)
	di `num_extract'
	
	forvalues i=1/`num_extract' { /* Loop through number of coefficients extracted */

		replace outcome = "`outcome'" if _n==`i'
		
		replace exp = "`exp'" if _n==`i'
		
		mat matobs = r(b)
		di "`: word `i' of `: colnames matobs''"
		replace coef = "`: word `i' of `: colnames matobs''" if _n==`i'
		
		replace nobs = r(N) if _n==`i'
		
		replace obs = matobs[1,`i'] if _n==`i'
		
		replace repsn = r(N_reps) if _n==`i'
		
		mat matct = r(c)
		replace cttrue = matct[1,`i'] if _n==`i'
		
		mat matreps = r(reps)
		replace repsnm = matreps[1,`i'] if _n==`i'
		
		mat matp = r(p)
		replace p = matp[1,`i'] if _n==`i'
		
		mat matci = r(ci)
		replace cilb = matci[1,`i'] if _n==`i'
		replace ciub = matci[2,`i'] if _n==`i'
		
		replace seed = `seed' if _n==`i'
		
		capture replace high_low = "`hl_name'" if _n==`i'
		capture replace dep = "`dep'" if _n==`i'
		
		capture replace base = `base' if _n==`i'
		
		}
		
		keep outcome-base
		
		keep if _n<=`num_extract'
	
*** Export data
	export delimited using "`outcome_csv'.csv", replace
	
	clear
	exit
