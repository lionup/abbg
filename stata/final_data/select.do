
clear 
clear matrix
clear mata
cap log close
set mem 200m
set more off
log using select, replace 

*** user parameters
global consumption  = 5  /*  5=totcons_psid_zero  *** The default                            
                             6=same as 5 without housing */
global min_wage     = 1  /*  1= Assigning missing for observations with wage lower then 1/2 state *** The default*/
                                                                                                                                                            
/* User parameter to control dropping outlyers */
 
global pec_drop     = 0.25  /* Drop observation in the bottom "pec_drop" percent of "jumps" */
                                                                                                                               
 
*********************** This sampling condition needs to be set also in the GMM file. Note that SEO are always dropped below.
*global sample2 "keep if age>=27 & age<=65"                                                          
                                                                                                                               
u data,clear
 
*** Liquid holdings
gen liquidhold = cash + bonds + stocks
 
* Adjust by price index
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
 
*${sample2}
 
gen log_ly  = log(ly) - log(price)
gen log_wly = log(wly) - log(price)

replace totcons = totcons/price 
gen log_totcons = log(totcons)	  /* log consumption */ 

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
  foreach var of varlist ly wly totcons {
    replace `var'=. if f2.pec_d_log_`var'_lag==1  /* Note that we are only assinging missing values to the year with the jump */
  }
  drop d_* pec_* log_*
}
 
gen totly = (ly+wly)/price    /* total family income */
gen log_totly  = log(totly) 
gen log_cons = log(totcons)	  /* log consumption */ 
 
gen wife_employed = (hourw>0 & hourw!=.)
gen mort1_dum = (mortgage1>0 & mortgage1!=.)
gen mort2_dum = (mortgage2>0 & mortgage2!=.)
 
log close
saveold select, replace 


