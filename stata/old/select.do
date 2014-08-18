
clear 
clear matrix
clear mata
cap log close
set mem 200m
set more off

*** user parameters
global consumption  = 5  /*  5=totcons_psid_zero  *** The default                            
                             6=same as 5 without housing */
global min_wage     = 1  /*  1= Assigning missing for observations with wage lower then 1/2 state *** The default*/
                                                                                                                                                            
/* User parameter to control dropping outlyers */
 
global pec_drop     = 0.25  /* Drop observation in the bottom "pec_drop" percent of "jumps" */
                                                                                                                               
 
*********************** This sampling condition needs to be set also in the GMM file. Note that SEO are always dropped below.
*global sample2 "keep if age>=27 & age<=65"                                                          
                                                                                                                               
*** see more parameters in Mata part
 
u data3,clear

cap log close
log using select, replace 
 
/*** Merge with welfare simulation data
ren year yearo
gen year = yearo+1900
sort id year
merge 1:1 id year using "$welfare/log_simulated_benefits"
drop if _m==2
drop _m
drop year
ren yearo year
 
egen   ben_temp = mean(all_benefits), by(year state kids)
replace all_benefits = ben_temp if all_benefits==.          /* Replace missing simulated benefits by state-year-kids mean */
*/
 
*** Liquid holdings
gen liquidhold = cash + bonds + stocks
 
*** merge prices
sort year
merge m:1 year using price_indices
keep if _m==3
drop _m
 
* Change base year to 2000
/*
foreach var of varlist p_* {
  su `var' if year==100
  local denom=r(mean)
  replace `var'=`var'/`denom'
}
*/
 
* Apply configuration
*gen price=p_totcons
replace tot_assets3 = tot_assets3/price
replace fin_assets1 = fin_assets1/price
 
sort person year
 
if $consumption == 6 {
  replace totcons = totcons_nh
}
 
* if hourly wage very small, measurement error
 if $min_wage==1 {
  gen oly = ly
  gen owly = wly
  replace ly=.  if (ly/hours)<0.5*min_wage
  replace wly=. if (wly/hourw)<0.5*min_wage
}
 
 
***************************
*** Additional Sampling ***
***************************
/* Baseline specification */
keep if seo==0
 
/* Want to use data from age 26 for lags at age 28, but don't want to condition the sample on characteristics at age 26 */
tsset person year
foreach var of varlist tot_assets* fin_assets* liquidhold {
  gen l2_`var'=l2.`var'
}
 
${sample2}
 
*gen no_work_d = (hours==0 | oly==0)  /*  Note that there are no missing hours. The missing ly is only if the minimum_wage criteria is applied, we assume that this is as in the case of ly==0 */
*gen no_work_d =(hours==0 | ly==0)   /*  Note that there are no missing hours. The missing ly is only if the minimum_wage criteria is applied, we assume that this is as in the case of ly==0 */
*egen todrop = max(no_work_d), by(person)
*drop if todrop == 1   /* note that we are only working in this setup with households where the prime earner is always working */
 
/* dependant variables and covariates for the estimation of the predicted part of wages, earnings and consumption as well as for selection equation */
 
*gen log_y  = log(ly) - log(price)
*gen log_yw = log(wly) - log(price)
*gen log_c  = log(totcons) - log(price)
*gen log_w  = log(ly/hours) - log(price)
*gen log_ww = log(wly/hourw) - log(price)
*gen log_toty= log(ly + wly) - log(price)                   /* This is only used to track inequality over time not in estimation */

gen totly = (ly+wly)/price    /* total family income */
gen log_totly  = log(totly)  /* log total family income */
 
replace totcons = totcons/price 
gen log_cons = log(totcons)	  /* log consumption */ 
 
gen wife_employed = (hourw>0 & hourw!=.)
gen mort1_dum = (mortgage1>0 & mortgage1!=.)
gen mort2_dum = (mortgage2>0 & mortgage2!=.)
 
if $pec_drop>0 {
 
  *** Use raw moments (not residual moments) that are associated with measurement error
  tsset person year
  local npec = 100/${pec_drop}
  di `npec'

  /* Generate the first difference for logs and the interaction between first difference and lagged first difference (which would be large in absolute value for large values of transitory shocks or for measurement error */
                 
  foreach var of varlist log_* {
    gen d_`var' = `var' - l2.`var'
    gen d_`var'_lag = d_`var'*l2.d_`var'
  }

  /* Generate percentiles of the interacted difference by year */
  foreach var of varlist d_*_lag {
    egen pec_`var'=  xtile(`var'), by(year) n(`npec')
  }

  /* Assign missing values for the variable with the potential measurement error  */
  foreach var of varlist log_*  {
    replace `var'=. if f2.pec_d_`var'_lag==1  /* Note that we are only assinging missing values to the year with the jump */
  }

  /* Check the new distribution of the interacted lag */
  foreach var of varlist log_* {
    gen d_`var'_trunc = `var' - l2.`var'
    gen d_`var'_trunc_lag = d_`var'_trunc*l2.d_`var'_trunc
  }
  su d_*_lag
  drop d_* pec_*

  *replace log_w = . if log_y==.
  *replace log_y = . if log_w==.

  *replace log_ww = . if log_yw==.
  *replace log_yw = . if log_ww==.
 
}
 
*drop if educ==. | weduc==.
log close
saveold select, replace 


