*** TAPBP VARIABLE CREATION

*** SETTINGS
    clear
	clear matrix
	clear mata
	set more off
	set maxvar 20000
    capture log close
	
 ** Create date local to name files
	local currdate_raw: display %td_CCYY_NN_DD date(c(current_date), "DMY")
	local currdate = subinstr(trim("`currdate_raw'"), " " , "-", .)
	
*** PATHS
	local path			= ""
	
	local datapath  	= "`path'/Data"
 	local logpath   	= "`path'/Log"
	local outputpath 	= "`path'/Output"

	local inputfile_KIR			= "feeney_malaria KIR"
	local inputfile_snp1		= "MalariaSNPAllTapasinDependence"
	local inputfile_snp2		= "MalariaVictoriafull032020[3]"
	local inputfile_HLA			= "#694 HLA typing results 7-18-17"
	local inputfile_EIR			= "CDC light trap data all 3 sites through June 2016 FINAL"
	local inputfile_ind			= "PRISM cohort study individual study participant database through June 30th 2016 FINAL"
	local inputfile_clin 		= "PRISM cohort study all clinic visits database through June 30th 2016 FINAL"
	local inputfile_langJIN		= "Jinja PRISM1 cohort - Primary languages.xls"
	local inputfile_langKAN		= "Kanungu PRISM1-primary languages.xlsx"
	local inputfile_langTOR 	= "Tororo PRISM1 cohort - Primary languages.xls"
	
	local outputfile	 		= "Malaria_TAPBP_Variables_`currdate'"
		
	local logfile	  			= "Malaria_TAPBP_Variables_`currdate'"

*** OPEN LOG FILE
	log using "`logpath'/`logfile'", replace
	
********************************************************************************

*** GET CLINICAL IDS FROM KIR DATASET
	import excel "`datapath'/`inputfile_KIR'.xlsx", sheet("Sheet1") firstrow case(lower) clear

	rename *, lower
	keep hgal pid

	drop if pid == "MF3273-4YNA"|pid == "MF3273-5A7P" /* UNSURE WHICH ARE THE TRUE RESULTS, DROP PER PERRI */
	drop if pid == "MF3287" /* NOT IN CLINICAL DATASET */
	
	replace pid = substr(pid, 3, 6)
	rename pid id
	destring id, replace

	save "`datapath'/`outputfile'", replace

********************************************************************************	

*** SNP DATA
	import excel "`datapath'/`inputfile_snp1'.xlsx", sheet("MalariaVictoria") firstrow case(lower) clear

	merge 1:1 hgal using "`datapath'/`outputfile'", gen(merge1)

	*drop if _m == 1 /* NO CLINICAL ID - HAVE TO DROP RESULT - 1 CASE */
	*drop if _m == 2 /* NO SNP RESULTS */
	
	*drop _m
	
	*keep hgal id rs* merge1 
	
*** Create variables
	gen snp_cg_gg = rs111686073=="CG"|rs111686073=="GG"
	tab snp_cg_gg rs111686073
	gen snp_ag_gg = rs59097151=="AG"|rs59097151=="GG"
	tab snp_ag_gg rs59097151
	
	gen snp_aa = rs59097151=="AA"
	tab snp_aa rs59097151
	
	save "`datapath'/`outputfile'", replace
	
*** SNP DATA: Tapasin
	import excel "`datapath'/`inputfile_snp2'.xlsx", sheet("MalariaVictoriafull032020") firstrow case(lower) clear

	merge 1:1 hgal using "`datapath'/`outputfile'", gen(merge2)

	drop if merge1==1|merge2==1 /* NO CLINICAL ID - HAVE TO DROP RESULT - 1 CASE */
	drop if merge1==2 & merge2==2 /* NO SNP RESULTS IN EITHER DATASET*/
	
	* Check IDs with only partially complete dependence data
	* browse if merge2==2 & (hlaadependence!="NA"|hlabdependence!="NA"|hlacdependence!="NA")
	
*** Create 3-level expression variable
	/* 	CC/GG and CC/AG (group 1, value to be assigned -0.5595136)
		CG/AG and CC/AA (group 2, value to be assigned 0.102363)
		CG/AA and GG/AA (group 3 value to be assigned 0.54373093) */
	
	* NOTE: confirmed that grouptapasinexpression and vicexpcat are equivalent
	* grouptapasinexpression has no missing values
	
	* Variable for everyone - even those missing dependence data
	gen exp_3_c_m = -0.5595136 if grouptapasinexpression==1
	replace exp_3_c_m = 0.102363 if grouptapasinexpression==2
	replace exp_3_c_m = 0.54373093 if grouptapasinexpression==3
	
	* Variable for those with complete dependence data
	gen exp_3_c = exp_3_c_m if vicexpcat<.
	
*** Create variables
	rename tapdep dep_all_o
	gen float dep_all = dep_all_o /* Make float variables for ritest */
	
	foreach x in a b c {
		rename `x'tap dep_hla_`x'_o
		gen float dep_hla_`x' = dep_hla_`x'_o /* Make float variables for ritest */
		}
	
	rename vicexpcat exp_3
	rename vicexpfullcat exp_6
	rename victapexp exp_t
	
*** Create variables: Missing dependence data, complete expression data
	gen exp_3_m = exp_3
		replace exp_3_m = grouptapasinexpression if exp_3_m==.
	gen exp_6_m = exp_6 
		replace exp_6_m = j if exp_6_m==.
	gen exp_t_m = exp_t
		replace exp_t_m = tapasinexpression if exp_t==.

*** Send SNP variables to missing if don't have complete dependence data
	gen snp_cg_gg_m = snp_cg_gg 
	gen snp_ag_gg_m = snp_ag_gg 
	gen snp_aa_m = snp_aa
	
	replace snp_cg_gg=. if merge2!=3 
	replace snp_ag_gg=. if merge2!=3 
	replace snp_aa=. if merge2!=3 

	keep id snp* exp* dep*
	tab1 snp* exp* 
	sum dep*
	
*** Histograms of dependence
/*	foreach x in all hla_a hla_b hla_c {
		hist dep_`x', name(hist_dep_`x', replace) width(0.1) start(0.3) freq addlabels xlabel(0.3(0.1)2.7, angle(90)) ylabel(0(50)250)
		graph export "`path'/Output_Graphs/hist_dep_`x'.pdf", replace
		}  */

*** Split dependence variables into high/low dichotomous
	/*  (Initial:
		Overall tapasin dependence: 1.9 
		HLA-A: 1.1
		HLA-B: 1.8
		HLA-C: 0.9)
		
		Final:
		1.70 for ABC, 1.10 for A, 1.40 for B and 0.80 for C*/
	
	local dep_split_all = 	"1.70"
	local dep_split_hla_a = "1.10"
	local dep_split_hla_b = "1.40"
	local dep_split_hla_c = "0.80"
	
	foreach x of local dep_split_hla_a {
		local varnum = substr("`i'", 1, 1) + substr("`i'", 3, 2)
		di "`varnum'"
		}
	
	label define highlow 0 "Low" 1 "High"
	
	foreach x in all hla_a hla_b hla_c {
		foreach i of local dep_split_`x' {
			* Create number with no decimal for variable name
			local varnum = substr("`i'", 1, 1) + substr("`i'", 3, 2)
			di "`varnum'"
			
			* Recode based on the list of splits
			recode dep_`x' (min/`i' = 0) (`i'/max = 1), gen (dep_`x'_`varnum'_b)
			label val dep_`x'_`varnum'_b highlow
			label var dep_`x'_`varnum'_b "dep_`x'_`varnum'_b"
			
			* Graph to check
			/*hist dep_`x', by(dep_`x'_`varnum'_b) name(hist_dep_`x'_`varnum'_b, replace) width(0.1) start(0.3) freq addlabels xlabel(0.3(0.1)2.7, angle(90)) ylabel(0(50)250)
			graph export "`path'/Output_Graphs/hist_dep_`x'_`varnum'_b.pdf", replace */
		}
	}

	
*** Create three level vars
	label define dep_3_a 0 "Independent" 1 "Dep+AG/GG" 2 "Dep+AA"
	
	foreach var of varlist dep_*_b {
		local x = regexr(regexr("`var'", "dep_", ""), "_b+$", "") 
			di "`x'"
		
		if "`x'" != "hla" { // Don't want the loops to use dep_hla_b
		gen dep_`x'_3_a = 0 if dep_`x'_b==0 /* dep LOW */
			replace dep_`x'_3_a = 1 if dep_`x'_b==1 & snp_ag_gg==1 /* AG/GG and dep HIGH */
			replace dep_`x'_3_a = 2 if dep_`x'_b==1 & snp_ag_gg==0  /* AA and dep HIGH */
		label val dep_`x'_3_a dep_3_a
		}
	
	}

	
	label define dep_3_c 0 "Independent" 1 "Dep+CC" 2 "Dep+CG/GG" // CG/GG = AA = high expression
	
	foreach var of varlist dep_*_b {
		local x = regexr(regexr("`var'", "dep_", ""), "_b+$", "") 
			di "`x'"
		
		if "`x'" != "hla" { // Don't want the loops to use dep_hla_b
		gen dep_`x'_3_c = 0 if dep_`x'_b==0 /* dep LOW */
			replace dep_`x'_3_c = 1 if dep_`x'_b==1 & snp_cg_gg==0 /* CC and dep HIGH */
			replace dep_`x'_3_c = 2 if dep_`x'_b==1 & snp_cg_gg==1  /* CG/GG and dep HIGH */
		label val dep_`x'_3_c dep_3_c
		}
	
	}

*** Order variables
	order *, alpha
	order id, first

*** Check sample sizes	
	foreach var of varlist * {
		qui count if `var'<.
		qui local count_n = r(N)
		di "`var'   `count_n'"
	}

	save "`datapath'/`outputfile'", replace

********************************************************************************
********************************************************************************

*** HLA VARIABLES

*** HLA
	import excel "`datapath'/`inputfile_HLA'.xlsx", sheet("Sheet1") firstrow case(lower) clear
	
	foreach var of varlist a_1-drb5_2 {
		gen `var's = substr(`var',1,2)
		destring `var's, replace
		
		gen hla_`var'_n = regexm(`var', "N") /* Tag if null allele */
		replace `var's=. if hla_`var'_n==1 /* Make missing if null allele */
		
		rename `var's hla_`var'_s
		gen `var'_length = length(`var')
		gen hla_`var'_l = `var'
		}

	drop if pid == "MF3273-4YNA"|pid == "MF3273-5A7P" /* UNSURE WHICH ARE THE TRUE RESULTS, DROP PER PERRI */
	drop if pid == "MF3287" /* NOT IN CLINICAL DATASET */

*** CREATE VARIABLES FOR HLA ALLELES
 ** FOUR DIGITS
	foreach h in hla_a hla_b hla_c {
		local hla_length_1 = substr("`h'", 5, .) + "_1_length"
		di "`hla_length_1'"
		
		local hla_length_2 = substr("`h'", 5, .) + "_2_length"
		di "`hla_length_2'"
		
		gen `h'_1_la=subinstr(`h'_1_l, ":", "", .) if `hla_length_1' > 2 /* remove : */
		gen `h'_2_la=subinstr(`h'_2_l, ":", "", .) if `hla_length_2' > 2 /* Send long version to missing if only 2 digits */
	
		levelsof `h'_1_la, local(levels) 			/* get levels of HLA */
		levelsof `h'_2_la, local(levels2)
		
		foreach l of local levels {
			gen `h'_l_`l'=.							/* create empty variable */
			}
		foreach l of local levels2 {
			capture confirm variable `h'_l_`l'
			if _rc {
				gen `h'_l_`l'=.						/* create var if not already present */
				local levels "`levels' `l'"			/* update local with levels missing from 1 */
			}
		}
		
		foreach l of local levels {
			* Make 1 if 4-digit allele 1 or 2 matches
			replace `h'_l_`l' = 1 if `h'_1_la=="`l'"|`h'_2_la=="`l'"
			* Make 0 if 4-digit allele 1 and 2 do not match and are not missing
			replace `h'_l_`l' = 0 if `h'_1_la!="`l'" & `h'_1_la!="" & `h'_2_la!="`l'" & `h'_2_la!=""
			* Make 0 if missing 4-digit allele, but non-missing 2-digit allele does not match
			replace `h'_l_`l' = 0 if `h'_l_`l'== . & ///
					(substr("`l'",1,2) != substr(`h'_1_l,1,2)|(`hla_length_1' > 2 & `h'_1_la!="`l'")) & `h'_1_l!= "" & ///
					(substr("`l'",1,2) != substr(`h'_2_l,1,2)|(`hla_length_2' > 2 & `h'_2_la!="`l'")) & `h'_2_l!= ""
			local hla_varlab = upper(substr("`h'",5,.)) + "*" + substr("`l'",1,2) + ":" + substr("`l'",3,.)
			label var `h'_l_`l' "`hla_varlab'"
			}
		browse `h'_*
	}
	
*** MERGE IN DATASET
	drop dpa1_1-hla_drb5_2_l
	replace pid = substr(pid, 3, 6)
	rename pid id
	destring id, replace
	merge m:1 id using "`datapath'/`outputfile'"
	rename _m merge_hla
	drop if merge_hla==1 /* Drop if only in HLA data */
	save "`datapath'/`outputfile'", replace

********************************************************************************
********************************************************************************

*** COHORT DATA

*** INDIVIDUAL - ONE RECORD PER INDIVIDUAL
	use "`datapath'/`inputfile_ind'", clear
	rename agecat agecat_enroll
	rename childadult child
	
 ** Male
	gen male=gender==1
	
 ** G6PD
	recode G6PD (0 1 3 = 0) (2 4 = 1) (9 = .), gen(g6pd_r)
	tab G6PD g6pd_r, m
	label var g6pd_r "G6PD: Male hemizygote/Female homozygote"
	
 ** HBS
	recode hbs (1 2= 1) (9 = .), gen(hbs_r)
	label var hbs_r "Hb AS/Hb SS"
	
 ** ALPHA THALASSEMIA
	recode alphathal (9 = .), gen(alphathal_r)
	label var alphathal_r "Alpha globin variant"
	label val alphathal_r alphathal

*** ALL STUDY VISITS - LONG FORM - MULTIPLE RECORDS PER INDIVIDUAL
	merge 1:m id using "`datapath'/`inputfile_clin'"
	drop _m
	
 ** DROP VISITS FOR TORORO (site 3) AFTER DEC 31, 2014
	drop if siteid==3 & date > date("20141231","YMD") /*CHECK DATE*/

*** INCIDENCE
 ** MALARIA INCIDENCE BY PERSON (ENTIRE STUDY)
  * PERSON-TIME FOR ENTIRE STUDY
	bysort id: egen firstvisit = min(date)
	bysort id: egen lastvisit = max(date)
	format %td firstvisit lastvisit
	gen pt_study = lastvisit-firstvisit
	browse id date firstvisit lastvisit pt_study
	
	gen pt_study_month = pt_study/(365/12)
	
  * EPISODES OF INCIDENT MALARIA
	bysort id: egen incidentmalaria_tot = total(incidentmalaria)
	
  * INCIDENCE RATE
	gen incidencerate = incidentmalaria_tot/pt_study
	
	browse id age date firstvisit lastvisit incidentmalaria* pt_study incidencerate
	
 ** MALARIA INCIDENCE BY PERSON BY YEAR
  * PERSON-TIME BY YEAR
   	* CREATE YEAR VARIABLE FOR EACH VISIT
	  gen year = year(date)
  
    * YEAR OF LAST VISIT
	  gen firstvisit_year = year(firstvisit)
	  
	* YEAR OF LAST VISIT
	  gen lastvisit_year = year(lastvisit)
	
	* PERSON-TIME PER YEAR
	  * Long form
	  gen pt_year = .
	  forvalues i = 2011/2016 {
		* for the first year a patient is in the study PT = year/12/31-first study visit
		replace pt_year = date("`i'1231","YMD")-firstvisit if firstvisit_year==`i' & year==`i'
		* for all in between years PT = 365 (except 2012 PT = 366: leap year)
		replace pt_year = 365 if firstvisit_year<`i' & lastvisit_year>`i' & year==`i'
		* for the last year a patient is in the study PT = last study visit - year/1/1
		replace pt_year = lastvisit-date("`i'0101","YMD") if lastvisit_year==`i' & year==`i'
		
		* add 1 to PT for 2012 due to leap year
		if `i' == 2012 { /* leap year */
			replace pt_year = pt_year + 1 if firstvisit_year<`i' & lastvisit_year>`i' & year==`i'
		}
	  }
	  
	  replace pt_year = . if firstvisit==lastvisit
	  replace pt_year = lastvisit-firstvisit if firstvisit_year==lastvisit_year
	  browse id date firstvisit lastvisit pt_*

	* EPISODES OF INCIDENT MALARIA BY YEAR
	  bysort id year: egen incidentmalaria_year = total(incidentmalaria)
	  browse id date incidentmalaria_year incidentmalaria
	
	* INCIDENCE RATE BY YEAR
	  gen incidencerate_year = incidentmalaria_year/pt_year
	  browse id age date firstvisit lastvisit year incidentmalaria_year pt_year incidencerate_year
	
	* CHECKS
	  /* keep id year pt_year incidencerate_year
		 duplicates drop
		 sum pt_year
		 sum incidencerate_year */

*** PARISITEMIA PREVALENCE (3 MONTH INTERVALS - BINARY Y/N)
	browse id date routinevisit age parasitedensity monthyear LAMP malariacat
	
 ** CREATE VISIT PERIOD
  * MIN DATE - 8/5/11
  * MAX DATE - 6/30/16
	gen month = month(date)

	gen q = 1 if month==1|month==2|month==3
		replace q = 2 if month==4|month==5|month==6
		replace q = 3 if month==7|month==8|month==9
		replace q = 4 if month==10|month==11|month==12
	
	bysort id year q: egen para_q_max = max(parasitedensity)
	bysort id year q: egen lamp_q = max(LAMP)
	
	recode para_q_max (0=0) (1/max=1), gen(para_q)
	
	gen para_lamp_q = 1 if (para_q==1|lamp_q==1)
		replace para_lamp_q = 0 if para_q==0 & (lamp_q==0|lamp_q==.)
		replace para_lamp_q = . if child==0
		/* Only look at this among children since LAMP not done on all adults */
		/* For now, there are 862 routine visits for which parasite density = 0 and LAMP is missing,
		   7% of visits that should have had a LAMP;
		   include these for now - this is a limitation of the analysis */	
	
	sort id date
	browse id visittype date month year q parasitedensity para_q*

*** LOG-TRANSFORMED PARASITEMIA
	gen para_ln = log(parasitedensity)

*** AGE 
 ** AGECAT
	rename agecat agecat_visit
	
 ** MEAN AGE FOR YEAR
	bysort id: egen age_first = min(age)
	bysort id: egen age_last = max(age)

	gen age_year_start = .
	gen age_year_end = .
	
	forvalues i = 2011/2016 {
		/* Calculate age as the midpoint of the year
           Midpoint of entire year if not first or last year observed
           Midpoint of time from first visit to Dec 31 for first year
           Midpoint of time from Jan 1 to last visit for last year */
		* for all years except first, make start age age on Jan 1
		replace age_year_start = (date("`i'0101","YMD")-dob)/365.25 if firstvisit_year!=`i' & year==`i'
		* for first year, make start age first recorded age
		replace age_year_start = age_first if firstvisit_year==`i' & year==`i'
		* for all years except last, make end age age on Dec 31
		replace age_year_end = (date("`i'1231","YMD")-dob)/365.25 if lastvisit_year!=`i' & year==`i'
		* for last year, make end age last recorded age
		replace age_year_end = age_last if lastvisit_year==`i' & year==`i'
		}
	
	gen age_y = (age_year_end+age_year_start)/2
	
	gen agecat_y = 1 if age_y<5
		replace agecat_y = 2 if age_y>=5 & age_y<12
		replace agecat_y = 3 if age_y>=18 & age_y<.
	tab age_y if agecat_y==.
	label val agecat_y agecat
	
 ** MEAN AGE PER QUARTER
	gen age_q_start = .
	gen age_q_end = .
	
	bysort id: gen firstvisit_q_raw = q if _n==1
	bysort id: egen firstvisit_q = max(firstvisit_q_raw)
	bysort id: gen lastvisit_q_raw = q if _n==_N
	bysort id: egen lastvisit_q = max(lastvisit_q_raw)
	
	forvalues i = 2011/2016 {
			/* 	Calculate age as the midpoint of the quarter
				Midpoint of entire quarter if not first or last quarter observed
				Midpoint of time from first visit to end of quarter for first quarter
				Midpoint of time from start of quarter to last visit for last quarter */
			replace age_q_start = (date("`i'0101","YMD")-dob)/365.25 if year==`i' & q==1
			replace age_q_start = (date("`i'0401","YMD")-dob)/365.25 if year==`i' & q==2
			replace age_q_start = (date("`i'0701","YMD")-dob)/365.25 if year==`i' & q==3
			replace age_q_start = (date("`i'1001","YMD")-dob)/365.25 if year==`i' & q==4
			replace age_q_start = age_first if firstvisit_year==`i' & year==`i' & firstvisit_q==q
			
			replace age_q_end = (date("`i'0331","YMD")-dob)/365.25 if year==`i' & q==1
			replace age_q_end = (date("`i'0630","YMD")-dob)/365.25 if year==`i' & q==2
			replace age_q_end = (date("`i'0930","YMD")-dob)/365.25 if year==`i' & q==3
			replace age_q_end = (date("`i'1231","YMD")-dob)/365.25 if year==`i' & q==4
			replace age_q_end = age_last if lastvisit_year==`i' & year==`i' & lastvisit_q==q
	}
	
	gen age_q = (age_q_end+age_q_start)/2
	
	gen agecat_q = 1 if age_q<5
		replace agecat_q = 2 if age_q>=5 & age_q<12
		replace agecat_q = 3 if age_q>=18 & age_q<.
	tab age_q if agecat_q==.
	label val agecat_q agecat
	
	browse id age year q age_q*
	
*** MERGE IN SNP DATASET
	merge m:1 id using "`datapath'/`outputfile'"
	rename _m merge_clinical
	
 ** MISSING FROM DATA		
	quietly tab id if merge_clinical==1
	display r(r)
	quietly tab id if merge_clinical==2
	display r(r)	

	keep if merge_clinical==3
	
	save "`datapath'/`outputfile'", replace

********************************************************************************

*** EIR: Household level
	use "`datapath'/`inputfile_EIR'", clear 

 ** Exclude time not in analysis
  * DROP VISITS FOR TORORO (site 3) AFTER DEC 31, 2014
	drop if siteid==3 & date > date("20141231","YMD") 

 ** Calculate P. falciparum Sporozoite rates by site
	bysort siteid: egen numberpositive_tot = total(numberpositive)
	bysort siteid: egen numbertested_tot = total(numbertested)
	gen pfsr = numberpositive_tot/numbertested_tot

 ** Calculate human biting rate: geometric average adding 0.5 to all values
	levelsof hhid, local(levels)
	
	gen hbr = .
	
	* Loop through all values of HHID /* Ameans doesn't work with egen */
	foreach l of local levels {
		ameans(totalanopheles) if hhid==`l', add(0.5)
		replace hbr = r(mean_g) if hhid==`l'
	}
	
 ** Calculate EIR: pfsr * hbr * 365 days
	gen eir = pfsr * hbr * 365
	
 ** Calculate log eir
	gen eir_ln = log(eir)
	
	keep hhid eir eir_ln
	duplicates drop
	
*** MERGE IN KIR/HLA DATASET
	merge 1:m hhid using "`datapath'/`outputfile'"
	rename _m merge_eir
	
	order id hhid siteid eir eir_ln

 ** MISSING FROM EIR DATA	
 	quietly tab id if merge_eir==3
	display r(r)
	*892 subjects are in EIR data, two IDs (1 HH) are missing EIR because it wasn't collected
	
	drop if hhid==201031206 /* Missing EIR */
	
	keep if merge_eir==3 /* EIR DATA COMPLETE */
	
	save "`datapath'/`outputfile'", replace
	
********************************************************************************

*** Language	
 ** Jinja
	import excel "`datapath'/`inputfile_langJIN'", sheet("Jinja cohort") firstrow case(lower) clear
	keep id language
	drop if id==1114 /* Inconsistent duplicate ID not in KIR data */
	
	tempfile languagedata    /* create a temporary file */
	save "`languagedata'"      /* save memory into the temporary file */
	
 ** Kanungu
	import excel "`datapath'/`inputfile_langKAN'", sheet("Kanungu") firstrow case(lower) clear
	keep id language
	
	append using "`languagedata'"
	save "`languagedata'", replace
	
 ** Kanungu
	import excel "`datapath'/`inputfile_langTOR'", sheet("Tororo cohort") firstrow case(lower) clear
	keep id language
	
	append using "`languagedata'"
	
	drop if id==3194 & language=="Muganda" /* Other HH members speak Jap */
			
 ** Merge to cohort data	
	merge 1:m id using "`datapath'/`outputfile'"
	keep if _m!=1 /* Drop if not in KIR data */
	save "`datapath'/`outputfile'", replace
	
 ** Create clean language variable
	gen language_r = lower(language)
	
  * Fill in using HH ID
	bysort hhid: egen languagemode = mode(language_r)
	replace language_r = languagemode if language_r==""
	drop languagemode
	
  * Clean language string
	replace language_r =  "luganda" if language_r == "ganda"
	replace language_r =  "ateso" if language_r == "ites0"
	replace language_r =  "ateso" if language_r == "iteso"
	replace language_r =  "ateso" if language_r == "itesot"
	replace language_r =  "dhopadhola " if language_r == "jaluo"
	replace language_r =  "dhopadhola " if language_r == "jap"
	replace language_r =  "rukiga" if language_r == "kiga"
	replace language_r =  "dhopadhola " if language_r == "ludama"
	replace language_r =  "rukiga" if language_r == "lukiga"
	replace language_r =  "runyankole" if language_r == "lunyankole"
	replace language_r =  "kinyarwanda" if language_r == "lunyarwanda"
	replace language_r =  "runyoro" if language_r == "lunyoro"
	replace language_r =  "samia" if language_r == "lusamya"
	replace language_r =  "ateso" if language_r == "luteso"
	replace language_r =  "lunyole" if language_r == "munyole"
	replace language_r =  "runyankole" if language_r == "nkore"
	replace language_r =  "kinyarwanda" if language_r == "rwandese"

 ** Ethnicity: 3 Category
	gen ethnicity = .
	replace ethnicity =  3 if language_r =="acholi"
	replace ethnicity =  2 if language_r =="ateso"
	replace ethnicity =  3 if language_r =="dhopadhola "
	replace ethnicity =  2 if language_r =="karamajong"
	replace ethnicity =  1 if language_r =="kikuyu"
	replace ethnicity =  1 if language_r =="kinyarwanda"
	replace ethnicity =  2 if language_r =="langi"
	replace ethnicity =  1 if language_r =="luganda"
	replace ethnicity =  3 if language_r =="lugbara"
	replace ethnicity =  1 if language_r =="lugisu"
	replace ethnicity =  1 if language_r =="lugwere"
	replace ethnicity =  1 if language_r =="lunyole"
	replace ethnicity =  1 if language_r =="lusoga"
	replace ethnicity =  1 if language_r =="lutoro"
	replace ethnicity =  1 if language_r =="rukiga"
	replace ethnicity =  1 if language_r =="runyankole"
	replace ethnicity =  1 if language_r =="runyoro"
	replace ethnicity =  1 if language_r =="samia"
	label define ethnicity 1 "Bantu" 2 "Nilo-Hamite" 3 "Luo"
	label val ethnicity ethnicity
	
 ** Bantu vs. non-Bantu
	recode ethnicity (2/3 = 0), gen(bantu)
	
 ** Site/Bantu
	gen site_bantu = siteid
		replace site_bantu = 0 if siteid==1 & bantu==1 
		replace site_bantu = . if siteid==1 & ethnicity==.
		label define site_bantu 0 "Jinja/Bantu" 1 "Jinja/Non-Bantu" 2 "Kanungu (99% Bantu)" 3 "Tororo (99% Non-Bantu)"
		label val site_bantu site_bantu

 ** Keep variables
	keep id hhid site* *bantu* snp* exp* dep* year para_ln incidentmalaria_year pt_year age age_y agecat_y q para_q age_q agecat_q male gender eir_ln g6pd_r hbs_r alphathal_r child year visittype parasitedensity febrile hla_b_l_5301 hla_c_l_0602
	drop dep_*_o

	save "`datapath'/`outputfile'", replace


