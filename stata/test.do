clear
use wen
drop _I*
xi i.educ i.yb i.weduc i.wyb
local varlist2 _I* bigcity 
************************************
* account for covariates effects in variances
* use method one with linear regression
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

sum uc muc utoty mutoty ua mua
drop uc2* utoty2* ua2*

* drop missing after standardization and construct balanced sample
drop if muc == . | mutoty==. | mua==.
by person, sort: gen numwav= _N
keep if numwav == 6
drop numwav
codebook person
