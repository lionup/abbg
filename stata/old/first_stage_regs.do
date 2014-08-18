
clear 
clear matrix
clear mata
cap log close
set more off
set mem 200m
log using first_stage_regs, replace

u select_2,clear
************************************
*** Balanced subsample N=749 T=6 ***
************************************
gen log_ass = log(tot_assets3)/* log asset */
*gen log_cons = log(totcons)	  /* log consumption */
*gen totly = ly+wly			  /* total family income */
*gen log_totly  = log(totly)  /* log total family income */

drop if log_ass == .		  /* drop if asset <= 0 or missing */
drop if log_cons == .         /* drop if consumption missing */
drop if log_totly == .		  /* drop if labor income missing */
drop if educ==. | weduc==.    /* drop if either male or female has education missing */
drop if state == .            /* drop if state missing */
drop if yb ==. | wyb == . 	  /* drop if either male or female has year of birth missing */
drop if kids == . 			  /* drop if # of kids missing */
drop if fsize ==.  			  /* drop if #of family member dummies missing */

* construct balanced sample
by person, sort: gen numwav= _N
keep if numwav == 6
drop numwav
codebook person

by year, sort: summarize totly totcons tot_assets3

*********************************
*** Regressions for residuals ***
*********************************
gen kidsout	= outkid==1       /* dummy for kids out*/
gen bigcity	= smsa==1|smsa==2 /* dummy for big city */
gen extra=(tyoth)>0           /* dummy for income recipient other than h/w */

/* Controls present:	
   ================
   year of birth dummies, education dummies, race dummies, # of kids dummies, 
   # of family member dummies, dummy for income recipient other than h/w, 
   state dummies, dummy for big city, dummy for kids out
   
   Interactions present:
   ====================
   1)male: educ * year of birth dummies 2)female: educ * year of birth dummies
*/	
xi	i.educ*i.yb i.weduc*i.wyb i.state i.fsize i.kids i.race i.wrace

* consumption regression 
reg log_cons kidsout bigcity extra _I*
predict uc if e(sample),res

* total income regression 
reg log_totly kidsout bigcity extra _I*
predict utoty if e(sample),res

* asset regression 
*reg tot_assets3 kidsout bigcity extra _I*
reg log_ass kidsout bigcity extra _I*
predict ua if e(sample),res

************************************
* account for covariates effects in variances
* use different set of controls
drop _I*
xi i.educ*i.yb i.weduc*i.wyb i.race i.wrace 
local varlist2 _I* bigcity 

* consumption
gen double uc2 = uc^2
reg uc2 `varlist2'
predict double uc2_yhat if e(sample)
gen muc = uc/sqrt(uc2_yhat)

* total income
gen double utoty2 = utoty^2
reg utoty2 `varlist2'
predict double utoty2_yhat if e(sample)
gen mutoty = utoty/sqrt(utoty2_yhat)

* total asset
gen double ua2 =ua^2
reg ua2 `varlist2'
predict double ua2_yhat if e(sample)
gen mua = ua/sqrt(ua2_yhat)

************************************
* if method one generate many negative values, 
* use nonlinear least square
*nl (utoty2 = {b0} * exp({xb: kidsout bigcity extra _I*}))
*predict double utoty2_yhat_nl if e(sample)
*gen mutoty_nl = utoty/sqrt(utoty2_yhat_nl)

*nl (ua2 = {b0} * exp({xb: kidsout bigcity extra _I*}))
*predict double ua2_yhat_nl if e(sample)
*gen mua_nl = ua/sqrt(ua2_yhat_nl)

sum uc muc utoty mutoty ua mua
drop uc2* utoty2* ua2*

* drop missing after standardization and construct balanced sample
*drop if muc == . | mutoty==. | mua==.
*by person, sort: gen numwav= _N
*keep if numwav == 6
*drop numwav
*codebook person

log close

saveold first_stage_resid_2, replace 


