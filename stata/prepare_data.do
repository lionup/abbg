
clear 
clear matrix
clear mata
cap log close
set mem 200m
set more off

*** user parameters
global participation = 0		/*  0=Shut off participation correction. 	*** The default
								1=with participation correction
								Note that we need to shut it off here since otherwise we estimate the residuals with the correction */
local  consumption 	=  5	/* 	1=nd + services								(this is the same as 5 since the sampling for 5 is already done in the imputation file).						
								2=nd 
								3=nd + services - childcare
								4=nd + services - rent 
								5=totcons_psid_zero							*** The default
								6=totcons_psid_miss
								7=totcons_cex_zero
								8=totcons_cex_miss
								9=totcons_zero_psid_onereg
								91=totcons_zero_psid_onereg2
								10=totcons_cex_onreg
								101=food
								102=fout
								103=gasoline
								104=electric
								105=water
								*/
global rhs=1				/* 	1=no hours or wages on RHS of imputation	*** The default: but doesn't matter since consumption is chosen to be 1 before
								2=only hours (no wages) on the RHS
								3=both hours and wages on the RHS
								*/ 
global assets		= 1		/* 	0= This period assets. 						*** The default
								1= Last period assets. */
global D2_part		= 1		/*  0= Calculate life time income of women as if they are always participating 
								1= Calculate life time income of women accounting for participation 				*** The default*/ 
global min_wage		= 1 	/* 	1= Assigning missing for observations with wage lower then 1/2 state minimum wage 	*** The default*/ 
global tau   		= 0		/* Linear tax rate */ 
					
					
/* User parameter to control dropping outlyers */ 

global pec_drop     = 0.25		/* Drop observation in the bottom "pec_drop" percent of "jumps" */ 
								

*** Baseline sample
global sample2 "keep if age>=30 & age<=65" 								
								
*** see more parameters in Mata part
di `consumption'
global consumption `consumption'
di $consumption

*** Tax schedule estimation
u totcons5_rhs1_part0,clear
do tax_rate_estimation

*********************************
*** Regressions for residuals ***
*********************************
tsset person year
tab year,gen(yrd)

/* head */ 
tab educ, gen(edd)
tab weduc, gen(wedd)
gen empl  	= empst==1
gen w_empl  = wempst==1
gen unempl	= empst==2
gen w_unempl= wempst==2
gen retir 	= empst==3
gen w_retir = wempst==3
gen white 	= race==1
gen w_white	= wrace==1 
gen black 	= race==2
gen w_black	= wrace==2 
gen other 	= race>=3
gen w_other = wrace>=3
gen selfe 	= self==3
gen w_selfe = selfw==3
tab state,gen(stated)
* gen veteran	=vet==1
* gen disability=disab==1
gen kidsout	= outkid==1
gen bigcity	= smsa==1|smsa==2
tab yb,gen(ybd)
tab wyb,gen(w_ybd)
tab fsize,gen(fd)
tab kids,gen(chd)
replace smsa=6 if smsa>6
* tab smsa,gen(cityd)
gen extra=(tyoth)>0
* tab age,gen(aged)
* tab agew,gen(w_aged)

cap log close
log using first_stage_regs, t replace

gen dlog_c=log_c-l2.log_c
gen dlog_w=log_w-l2.log_w
gen dlog_ww=log_ww-l2.log_ww
gen dlog_y=log_y-l2.log_y
gen dlog_yw=log_yw-l2.log_yw
gen dlog_toty=log_toty-l2.log_toty

foreach var of varlist kidsout bigcity kids fsize empl w_empl unempl w_unempl retir w_retir extra {
	gen d`var'=`var'-l2.`var'
}

tab dkids, gen(dchd) 
tab dfsize, gen(dfd)

/* consumption regression */ 

/* Controls present:	
   ================
   year dummies, year of birth dummies, education dummies, race dummies,# of kids dummies (level and change), 
   #of family member dummies (level and change), empl. status dummies at time of interview (level and change), 
   dummy for income recipient other than h/w (level and change), state dummies, dummy for big city (level and change), 
   dummy for kids not in FU (level and change) 
   
   Interactions present:
   ====================
   1)educ*year 2)race dummies*year 3)empl. status dummies*year (level and change) 4)state*year EXCLUDED 
   5)big city*year (level and change) 
   
*/	

#delimit;
xi	i.educ*i.year i.weduc*i.year i.white*i.year i.w_white*i.year i.black*i.year i.w_black*i.year 
	i.other*i.year i.w_other*i.year i.bigcity*i.year i.dbigcity*i.year 
	i.empl*i.year i.unempl*i.year i.retir*i.year i.w_empl*i.year 
	i.w_unempl*i.year i.w_retir*i.year i.dempl*i.year i.dunempl*i.year i.dretir*i.year 
	i.dw_empl*i.year i.dw_unempl*i.year i.dw_retir*i.year
	i.state*i.year i.mort1_dum*i.year i.mort2_dum*i.year;
reg dlog_c   	yrd* ybd* w_ybd* edd* wedd* white w_white black w_black other w_other 
				fd* chd* empl  unempl retir w_empl  w_unempl w_retir kidsout bigcity extra 
				dfd* dchd* dempl dunempl dretir dw_empl dw_unempl dw_retir dextra dkidsout dbigcity dextra 
				stated* 
				_IeduXyea_* _IwedXyea_* _IwhiXyea_* _Iw_wXyea_* _IblaXyea_* _Iw_bXyea_*
				_IothXyea_* _Iw_other_* 
				_IempXyea_* _Iw_empl_* _IuneXyea_* _Iw_unempl_* _IretXyea_* _Iw_retir_*
				_IbigXyea_* _IdbiXyea_*;
predict uc if e(sample),res;
#delimit cr

/* Earnings regression head - same covariates as consumption */ 
#delimit;
reg dlog_y   	yrd* ybd* w_ybd* edd* wedd* white w_white black w_black other w_other 
				fd* chd* empl  unempl retir w_empl  w_unempl w_retir kidsout bigcity extra 
				dfd* dchd* dempl dunempl dretir dw_empl dw_unempl dw_retir dextra dkidsout dbigcity dextra 
				stated* 
				_IeduXyea_* _IwedXyea_* _IwhiXyea_* _Iw_wXyea_* _IblaXyea_* _Iw_bXyea_*
				_IothXyea_* _Iw_other_* 
				_IempXyea_* _Iw_empl_* _IuneXyea_* _Iw_unempl_* _IretXyea_* _Iw_retir_*
				_IbigXyea_* _IdbiXyea_*;
predict uy if e(sample),res;
#delimit cr

/* Earnings regression Total household earnings - THIS IS ONLY USED TO FOLLOW EVOLUTION OF CONSUMPTION IN EQUALITY */ 
#delimit;
reg dlog_toty   yrd* ybd* w_ybd* edd* wedd* white w_white black w_black other w_other 
				fd* chd* empl  unempl retir w_empl  w_unempl w_retir kidsout bigcity extra 
				dfd* dchd* dempl dunempl dretir dw_empl dw_unempl dw_retir dextra dkidsout dbigcity dextra 
				stated* 
				_IeduXyea_* _IwedXyea_* _IwhiXyea_* _Iw_wXyea_* _IblaXyea_* _Iw_bXyea_*
				_IothXyea_* _Iw_other_* 
				_IempXyea_* _Iw_empl_* _IuneXyea_* _Iw_unempl_* _IretXyea_* _Iw_retir_*
				_IbigXyea_* _IdbiXyea_*;
predict utoty if e(sample),res;
#delimit cr

/* Wage regression head  */ 

/* Controls present:	
   ================
   year dummies, year of birth dummies, education dummies, race dummies, state dummies, 
   dummy for big city (level and change)
   
   Interactions present: 
   ====================
   1)educ*year 2)race dummies*year 3)state*year EXCLUDED 4)big city*year 

*/	
   
#delimit;
reg dlog_w	    yrd* ybd* edd*  white black other stated* bigcity dbigcity
			    _IeduXyea_* _IwhiXyea_* _IblaXyea_* _IothXyea_* _IbigXyea_* _IdbiXyea_*;
predict uw if e(sample),res;
#delimit cr

/* selection equation for the wife - same controls as consumption plus:
   1) stateXyear interaction 
   2) mortgageXyear interaction for first and second mortgage
*/   
   
#delimit;
probit wife_employed   	yrd* ybd* w_ybd* edd* wedd* white w_white black w_black other w_other 
				fd* chd* empl  unempl retir w_empl  w_unempl w_retir kidsout bigcity extra 
				dfd* dchd* dempl dunempl dretir dw_empl dw_unempl dw_retir dextra dkidsout dbigcity dextra 
				stated* 
				_IeduXyea_* _IwedXyea_* _IwhiXyea_* _Iw_wXyea_* _IblaXyea_* _Iw_bXyea_*
				_IothXyea_* _Iw_other_* 
				_IempXyea_* _Iw_empl_* _IuneXyea_* _Iw_unempl_* _IretXyea_* _Iw_retir_*
				_IbigXyea_* _IdbiXyea_*
				_IstaXyea_* _Imor*;
#delimit cr
testparm _IstaXyea_*
testparm _Imor*
testparm _IstaXyea_* _Imor*

predict z_eta if e(sample), xb
gen pdf = normalden(z_eta)
gen cdf = normal(z_eta)
gen inv_mills = pdf/cdf
drop pdf cdf

if ${participation}==0 {
	replace inv_mills=0		/* This would imply that it is ommited from the regression in what follows */ 
}

/* Earnings regression wife - same controls as consumption plus inverse mills ratio if selection correction applied */ 
gen dinv_mills=inv_mills-l2.inv_mills

#delimit;
reg dlog_yw   	yrd* ybd* w_ybd* edd* wedd* white w_white black w_black other w_other 
				fd* chd* empl  unempl retir w_empl  w_unempl w_retir kidsout bigcity extra 
				dfd* dchd* dempl dunempl dretir dw_empl dw_unempl dw_retir dextra dkidsout dbigcity dextra 
				stated* 
				_IeduXyea_* _IwedXyea_* _IwhiXyea_* _Iw_wXyea_* _IblaXyea_* _Iw_bXyea_*
				_IothXyea_* _Iw_other_* 
				_IempXyea_* _Iw_empl_* _IuneXyea_* _Iw_unempl_* _IretXyea_* _Iw_retir_*
				_IbigXyea_* _IdbiXyea_*
				dinv_mills;
local gamma=_b[dinv_mills];
predict uyw_temp if e(sample),res;
gen uyw = uyw_temp + `gamma'*dinv_mills;
drop uyw_temp;
#delimit cr

/* Wage regression wife  */ 

/* Controls present:
   ================
   year dummies, year of birth dummies, education dummies, race dummies, state dummies,
   dummy for big city (level and change)

   Interactions present: 
   ====================
   1)educ*year 2)race dummies*year 3)state*year EXCLUDED 4)big city*year 
   
*/	
   
#delimit;
reg dlog_ww 	yrd* w_ybd* wedd*  w_white w_black w_other stated* bigcity dbigcity
				_IwedXyea_* _Iw_wXyea_* _Iw_bXyea_* _Iw_oXyea_* _IbigXyea_* _IdbiXyea_*   
				dinv_mills;
local gamma=_b[dinv_mills];
predict uww_temp if e(sample),res;
gen uww = uww_temp + `gamma'*dinv_mills;
drop uww_temp;
#delimit cr

drop _I* 
log close

save data3_resid, replace 

/* Generate variables with explicit expressions for pi and for s */

/*  ************************* NOTEs ****************************
	1) 	Need to fix T. The sample only goes to 65, so for now I
		put 66 at the end of the sum, so if you start at 65 you 
		do not have zero (otherwise I might get zeros in the 
		denominator). 
	2) 	Do year by year or with year dummies? Year by year is 
		easier for the code below.
	3) 	Especially if doing it yearly, need to turn into a 
		function.
	***********************************************************/ 
	
*** User Parameters ***
scalar r	= 0.02
scalar tau	= $tau
scalar T 	= 65

u data3_resid, clear
gen log_h = log(hours)		/* log hours is only going to be used for calculation of ME variance */
gen log_hw = log(hourw)


/* The following calculates D1 */ 

gen age2 = age^2
gen age3 = age^3
xi	i.educ*age i.educ*age2 i.educ*age3 i.race*age i.race*age2 i.race*age3 i.yb 

reg log_y  age age2 age3 _I*
*predict test
*gen exp_test=exp(test) /* for debug - can chage check that the same as D1 for age 65 */ 

matrix beta_male = e(b) 

gen D1 = .

su yb
local yb_min = r(min)+1	/* Since first dummy is ommited */ 
local yb_max = r(max)
local jj_min = `yb_min'+1

mata
r		= st_numscalar("r")
tau		= st_numscalar("tau")
T		= st_numscalar("T")

yb_string = " _Iyb_" + st_local("yb_min")
jj = strtoreal(st_local("jj_min"))

while (jj<=strtoreal(st_local("yb_max"))) {
	yb_string = yb_string + " _Iyb_" + strofreal(jj)
	jj = jj + 1
}

data 	= st_data(.,("age","age2","age3","_Ieduc_2","_Ieduc_3", ///
					"_IeduXage_2","_IeduXage_3","_IeduXage2_2","_IeduXage2_3","_IeduXage3_2", ///
					"_IeduXage3_3","_Irace_2","_Irace_3","_IracXage_2","_IracXage_3","_IracXage2_2", ///
					"_IracXage2_3","_IracXage3_2","_IracXage3_3",yb_string))
				
ones 	= J(rows(data),1,1)
data 	= data,ones
beta_male 	= st_matrix("beta_male")
D1		= J(rows(data),1,0)

data1	= data	/* These are used later for the non-linear tax case */
beta1	= beta_male

i=1
j=0
while (i>0) {
	index	= data[.,1]:<T+1				/* assigns 1 for observations which are less or equal to 65 in the loop */ 
	i		= max(index)
	D1 		= D1 + exp((data*beta_male') - j*log(1+r)*ones):*index
	data[.,1]	= data[.,1] + ones
	data[.,2]	= data[.,1]:^2
	data[.,3]	= data[.,1]:^3
	data[.,6]	= data[.,1]:*data[.,4]
	data[.,7]	= data[.,1]:*data[.,5]
	data[.,8]	= data[.,2]:*data[.,4]
	data[.,9]	= data[.,2]:*data[.,5]
	data[.,10]	= data[.,3]:*data[.,4]
	data[.,11]	= data[.,3]:*data[.,5]
	data[.,14]	= data[.,1]:*data[.,12]
	data[.,15]	= data[.,1]:*data[.,13]
	data[.,16]	= data[.,2]:*data[.,12]
	data[.,17]	= data[.,2]:*data[.,13]
	data[.,18]	= data[.,3]:*data[.,12]
	data[.,19]	= data[.,3]:*data[.,13]
											/* Note that no need to update anything that does not have age in it */ 
	j = j+1
}

D1 = (1-tau)*D1
st_store(.,"D1",D1)
end
drop _I*

/* The following calculates D1 in an alternative way, using regressions on growth not levels */ 
gen l2_log_y=l2.log_y

xi	i.educ*age i.educ*age2 i.educ*age3 i.race*age i.race*age2 i.race*age3 i.yb 
drop _Iyb_33 _Iyb_34 _Iyb_76 _Iyb_77

reg dlog_y  age age2 age3 _I*	/* Note: a regression of dlog_y not log_y */

*predict test
*gen exp_test=exp(test) /* for debug - can chage check that the same as D1 for age 65 */ 

matrix beta_male = e(b) 

gen D1_alt = .

su yb if dlog_y!=.
local yb_min = r(min)+1	/* Since first dummy is ommited */ 
local yb_max = r(max)
local jj_min = `yb_min'+1

mata
r		= st_numscalar("r")
tau		= st_numscalar("tau")
T		= st_numscalar("T")

yb_string = " _Iyb_" + st_local("yb_min")
jj = strtoreal(st_local("jj_min"))

while (jj<=strtoreal(st_local("yb_max"))) {
	yb_string = yb_string + " _Iyb_" + strofreal(jj)
	jj = jj + 1
}

log_y	= st_data(.,"l2_log_y")
data 	= st_data(.,("age","age2","age3","_Ieduc_2","_Ieduc_3", ///
					"_IeduXage_2","_IeduXage_3","_IeduXage2_2","_IeduXage2_3","_IeduXage3_2", ///
					"_IeduXage3_3","_Irace_2","_Irace_3","_IracXage_2","_IracXage_3","_IracXage2_2", ///
					"_IracXage2_3","_IracXage3_2","_IracXage3_3",yb_string))
				
ones 	= J(rows(data),1,1)
data 	= data,ones
beta_male 	= st_matrix("beta_male")
D1_alt		= J(rows(data),1,0)

i=1
j=0
while (i>0) {
	index	= data[.,1]:<T+1				/* assigns 1 for observations which are less or equal to 65 in the loop */ 
	i		= max(index)
	log_y   = log_y + (data*beta_male')		/* Projecting log_y using log_y at t-1 */ 
	D1_alt 		= D1_alt + exp(log_y - j*log(1+r)*ones):*index
	data[.,1]	= data[.,1] + ones
	data[.,2]	= data[.,1]:^2
	data[.,3]	= data[.,1]:^3
	data[.,6]	= data[.,1]:*data[.,4]
	data[.,7]	= data[.,1]:*data[.,5]
	data[.,8]	= data[.,2]:*data[.,4]
	data[.,9]	= data[.,2]:*data[.,5]
	data[.,10]	= data[.,3]:*data[.,4]
	data[.,11]	= data[.,3]:*data[.,5]
	data[.,14]	= data[.,1]:*data[.,12]
	data[.,15]	= data[.,1]:*data[.,13]
	data[.,16]	= data[.,2]:*data[.,12]
	data[.,17]	= data[.,2]:*data[.,13]
	data[.,18]	= data[.,3]:*data[.,12]
	data[.,19]	= data[.,3]:*data[.,13]
											/* Note that no need to update anything that does not have age in it */ 
	j = j+1
}

D1_alt = (1-tau)*D1_alt
st_store(.,"D1_alt",D1_alt)
end
drop _I*


/* The following piece of code calculates D2 */ 

/* 	*****************************************************************************
	Notes: 
	1) 	The Mills ratio from the probit above is used in the wage equations. 
		These estimates are NOT going to be used in the prediction of the probability 
		for participation at different ages since it has covariates that change over
		time. However, I need to use this equation to get a consistent estimate for 
		the Mills ratio (it has the exclusion restrictions). 
	2) 	The earnings equation is calculated the same as for male, only using 
		also the Mills ratio. The prediction will NOT use the mills ratio (i.e.
		it is kept as part of the residual as before)
	3) 	Unlike D1, D2 might be equal zero (since we do not sample by wife's age).
		This would happen if the wife is older than the husband and is older than 
		65 at time t. 
	*****************************************************************************/

/* Estimate wage equations for wife */
replace wrace = 3 if wrace>3
gen agew2 = agew^2
gen agew3 = agew^3

*Note in this case we take 5 years birth cohorts to avoid ommitted cohorts 
gen wyb_coh=1
replace wyb_coh=2 if wyb>35 	/* First cell is 10 years since it is very small */ 
replace wyb_coh=3 if wyb>40
replace wyb_coh=4 if wyb>45
replace wyb_coh=5 if wyb>50
replace wyb_coh=6 if wyb>55
replace wyb_coh=7 if wyb>60
replace wyb_coh=8 if wyb>65
replace wyb_coh=9 if wyb>70
replace wyb_coh=10 if wyb>75 	/* Last cell is 10 years since it is very small */ 

xi	i.weduc*agew i.weduc*agew2 i.weduc*agew3 i.wrace*agew i.wrace*agew2 i.wrace*agew3 i.wyb_coh

reg log_yw  agew agew2 agew3 _I* inv_mills

/*

ren inv_mills inv_mills_old
gen inv_mills=0
predict test
drop inv_mills
ren inv_mills_old inv_mills
gen exp_test=exp(test)

*/

matrix beta_female = e(b) 

/* Estimate participation equations for wife using only age and fixed covariates */ 
probit wife_employed  agew agew2 agew3 _I*
matrix beta_prob = e(b) 

gen D2 = .

*su wyb		/* If cohorts are by year */
su wyb_coh
local yb_max = r(max)

mata
r		= st_numscalar("r")
tau		= st_numscalar("tau")
T		= st_numscalar("T")

// yb_string = " _Iwyb_26 _Iwyb_29 _Iwyb_31"	// If cohorts are by year have missing variable. Fix before use
// jj = 33

yb_string = ""
jj = 2

while (jj<=strtoreal(st_local("yb_max"))) {
	yb_string = yb_string + " _Iwyb_coh_" + strofreal(jj)
	jj = jj + 1
}

data 	= st_data(.,("agew","agew2","agew3","_Iweduc_2","_Iweduc_3", ///
					"_IwedXagew_2","_IwedXagew_3","_IwedXagew2_2","_IwedXagew2_3","_IwedXagew3_2", ///
					"_IwedXagew3_3","_Iwrace_2","_Iwrace_3","_IwraXagew_2","_IwraXagew_3","_IwraXagew2_2", ///
					"_IwraXagew2_3","_IwraXagew3_2","_IwraXagew3_3",yb_string))
						
ones 	= J(rows(data),1,1)
data 	= data,ones

beta_female = st_matrix("beta_female")
beta_female	= beta_female[1..(cols(beta_female)-2)],beta_female[cols(beta_female)] // Since we don't want the coefficient on the Mills ratio
beta_prob = st_matrix("beta_prob")

data2	= data	/* These are used later for the non-linear tax case */
beta2	= beta_female
beta2p	= beta_prob

D2		= J(rows(data),1,0)

D2_part = st_global("D2_part")

i=1
j=0
while (i>0) {
	index	= data[.,1]:<T+1
	i		= max(index)
	P		= normal(data*beta_prob')
	if (D2_part!= "1") {
			P=P:^0		/* This is the case where we do not account for participation when calculating expected life time income */ 
	}
	D2 		= D2 + exp(log(exp(data*beta_female'):*P) - j*log(1+r)*ones):*index
	data[.,1]	= data[.,1] + ones
	data[.,2]	= data[.,1]:^2
	data[.,3]	= data[.,1]:^3
	data[.,6]	= data[.,1]:*data[.,4]
	data[.,7]	= data[.,1]:*data[.,5]
	data[.,8]	= data[.,2]:*data[.,4]
	data[.,9]	= data[.,2]:*data[.,5]
	data[.,10]	= data[.,3]:*data[.,4]
	data[.,11]	= data[.,3]:*data[.,5]
	data[.,14]	= data[.,1]:*data[.,12]
	data[.,15]	= data[.,1]:*data[.,13]
	data[.,16]	= data[.,2]:*data[.,12]
	data[.,17]	= data[.,2]:*data[.,13]
	data[.,18]	= data[.,3]:*data[.,12]
	data[.,19]	= data[.,3]:*data[.,13]

	j = j+1
}

D2 = (1-tau)*D2
st_store(.,"D2",D2)

end


/* An alternative calculation for D2 that uses earnings at t-2 */ 

/* Estimate wage equations for wife */

xi	i.weduc*agew i.weduc*agew2 i.weduc*agew3 i.wrace*agew i.wrace*agew2 i.wrace*agew3 i.wyb_coh

gen l2_log_yw=l2.log_yw

*** If want to take out the cohort variables 
foreach var of varlist  _Iwyb_coh_2- _Iwyb_coh_10 {
	replace `var'=0 if `var'!=.
}	

reg dlog_yw  agew agew2 agew3 _I* dinv_mills

matrix beta_female = e(b) 

/* Estimate participation equations for wife using only age and fixed covariates */ 
probit wife_employed  agew agew2 agew3 _I*
matrix beta_prob = e(b) 

gen D2_alt = .

*su wyb		/* If cohorts are by year */
su wyb_coh
local yb_max = r(max)

mata
r		= st_numscalar("r")
tau		= st_numscalar("tau")
T		= st_numscalar("T")

// yb_string = " _Iwyb_26 _Iwyb_29 _Iwyb_31"	// If cohorts are by year have missing variable. Fix before use
// jj = 33

yb_string = ""
jj = 2

while (jj<=strtoreal(st_local("yb_max"))) {
	yb_string = yb_string + " _Iwyb_coh_" + strofreal(jj)
	jj = jj + 1
}

log_yw	= st_data(.,"l2_log_yw")
data 	= st_data(.,("agew","agew2","agew3","_Iweduc_2","_Iweduc_3", ///
					"_IwedXagew_2","_IwedXagew_3","_IwedXagew2_2","_IwedXagew2_3","_IwedXagew3_2", ///
					"_IwedXagew3_3","_Iwrace_2","_Iwrace_3","_IwraXagew_2","_IwraXagew_3","_IwraXagew2_2", ///
					"_IwraXagew2_3","_IwraXagew3_2","_IwraXagew3_3",yb_string))
ones 	= J(rows(data),1,1)
data 	= data,ones
beta_female = st_matrix("beta_female")
beta_female	= beta_female[1..(cols(beta_female)-2)],beta_female[cols(beta_female)] // Since we don't want the coefficient on the Mills ratio
beta_prob = st_matrix("beta_prob")
D2_alt		= J(rows(data),1,0)

D2_part = st_global("D2_part")

i=1
j=0
while (i>0) {
	index	= data[.,1]:<T+1
	i		= max(index)
	P		= normal(data*beta_prob')
	if (D2_part!= "1") {
			P=P:^0		/* This is the case where we do not account for participation when calculating expected life time income */ 
	}
	log_yw   = log_yw + (data*beta_female')		/* Projecting log_yw using log_yw at t-1. Note that this makes pi conditioned on participation */ 
	D2_alt 		= D2_alt + exp(log(exp(log_yw):*P) - j*log(1+r)*ones):*index
	data[.,1]	= data[.,1] + ones
	data[.,2]	= data[.,1]:^2
	data[.,3]	= data[.,1]:^3
	data[.,6]	= data[.,1]:*data[.,4]
	data[.,7]	= data[.,1]:*data[.,5]
	data[.,8]	= data[.,2]:*data[.,4]
	data[.,9]	= data[.,2]:*data[.,5]
	data[.,10]	= data[.,3]:*data[.,4]
	data[.,11]	= data[.,3]:*data[.,5]
	data[.,14]	= data[.,1]:*data[.,12]
	data[.,15]	= data[.,1]:*data[.,13]
	data[.,16]	= data[.,2]:*data[.,12]
	data[.,17]	= data[.,2]:*data[.,13]
	data[.,18]	= data[.,3]:*data[.,12]
	data[.,19]	= data[.,3]:*data[.,13]

	j = j+1
}

D2_alt = (1-tau)*D2_alt
st_store(.,"D2_alt",D2_alt)

end


/* The following calcualtes Q1 (this is the non-linear tax case) */
/* Note: It is an empirical question how to project total income. A robustness test would be 
   we use the projection for each earner separately instead of the projection of the sum. */ 

xi	i.educ*age i.educ*age2 i.educ*age3 i.race*age i.race*age2 i.race*age3 i.yb i.weduc*agew i.weduc*agew2 i.weduc*agew3 i.wrace*agew i.wrace*agew2 i.wrace*agew3 i.wyb_coh
gen log_toty_at = (1-tau)*log_toty + log(xsi)

reg log_toty_at age age2 age3 agew agew2 agew3 _I* 
*predict test
*gen exp_test=exp(test) /* for debug - can chage check that the same as D1 for age 65 */ 

matrix beta_tot = e(b) 

gen Q1 = .

su yb
local yb_min = r(min)+1	/* Since first dummy is ommited */ 
local yb_max = r(max)
local jj_min = `yb_min'+1

su wyb_coh
local wyb_max = r(max)

mata
r		= st_numscalar("r")
tau		= st_numscalar("tau")
T		= st_numscalar("T")

yb_string = " _Iyb_" + st_local("yb_min")
jj = strtoreal(st_local("jj_min"))

while (jj<=strtoreal(st_local("yb_max"))) {
	yb_string = yb_string + " _Iyb_" + strofreal(jj)
	jj = jj + 1
}

wyb_string = ""
jj = 2

while (jj<=strtoreal(st_local("wyb_max"))) {
	wyb_string = wyb_string + " _Iwyb_coh_" + strofreal(jj)
	jj = jj + 1
}


data 	= st_data(.,("age","age2","age3", "agew","agew2","agew3", /// 
					"_Ieduc_2","_Ieduc_3", "_IeduXage_2","_IeduXage_3","_IeduXage2_2","_IeduXage2_3","_IeduXage3_2", ///
					"_IeduXage3_3","_Irace_2","_Irace_3","_IracXage_2","_IracXage_3","_IracXage2_2", ///
					"_IracXage2_3","_IracXage3_2","_IracXage3_3",yb_string, /// 
					"_Iweduc_2","_Iweduc_3", "_IwedXagew_2","_IwedXagew_3","_IwedXagew2_2","_IwedXagew2_3", ///
					"_IwedXagew3_2", "_IwedXagew3_3","_Iwrace_2","_Iwrace_3","_IwraXagew_2","_IwraXagew_3","_IwraXagew2_2", ///
					"_IwraXagew2_3","_IwraXagew3_2","_IwraXagew3_3",wyb_string))
				
ones 	= J(rows(data),1,1)
data 	= data,ones

beta_tot 	= st_matrix("beta_tot")
Q1		= J(rows(data),1,0)

data3	= data
beta3	=beta_tot

i=1
j=0

/* index assigns 1 for observations for which the head is less or equal to 65 in the loop */

while (i>0) {
	index	= data[.,1]:<T+1				 
	i		= max(index)
	Q1 		= Q1 + exp((data*beta_tot') - j*log(1+r)*ones):*index
	data[.,1]	= data[.,1] + ones
	data[.,2]	= data[.,1]:^2
	data[.,3]	= data[.,1]:^3
	data[.,4]	= data[.,4] + ones
	data[.,5]	= data[.,4]:^2
	data[.,6]	= data[.,4]:^3
	data[.,9]	= data[.,1]:*data[.,7]
	data[.,10]	= data[.,1]:*data[.,8]
	data[.,11]	= data[.,2]:*data[.,7]
	data[.,12]	= data[.,2]:*data[.,8]
	data[.,13]	= data[.,3]:*data[.,7]
	data[.,14]	= data[.,3]:*data[.,8]
	data[.,17]	= data[.,1]:*data[.,15]
	data[.,18]	= data[.,1]:*data[.,16]
	data[.,19]	= data[.,2]:*data[.,15]
	data[.,20]	= data[.,2]:*data[.,16]
	data[.,21]	= data[.,3]:*data[.,15]
	data[.,22]	= data[.,3]:*data[.,16]
	data[.,70]	= data[.,4]:*data[.,68]
	data[.,71]	= data[.,4]:*data[.,69]
	data[.,72]	= data[.,5]:*data[.,68]
	data[.,73]	= data[.,5]:*data[.,69]
	data[.,74]	= data[.,6]:*data[.,68]
	data[.,75]	= data[.,6]:*data[.,69]
	data[.,78]	= data[.,4]:*data[.,76]
	data[.,79]	= data[.,4]:*data[.,77]
	data[.,80]	= data[.,5]:*data[.,76]
	data[.,81]	= data[.,5]:*data[.,77]
	data[.,82]	= data[.,6]:*data[.,76]
	data[.,83]	= data[.,6]:*data[.,77]
											/* Note that no need to update anything that does not have age in it */ 
	j = j+1
}

st_store(.,"Q1",Q1)
end
drop _I*











/* The following calcualtes s_hat (this is the non-linear tax case) */
/* Note: This is calculated as the sum of alpha*q_tilde rather than assuming that alpha
   is not correlated with q_tilde */ 

gen s_hat = .

mata
s_hat	= J(rows(data1),1,0)
i=1
j=0

* Note that the timing has to match the timing in the equation - q it at t-1
l1_data1=data1
l1_data1[.,1]	= data1[.,1] - ones
l1_data1[.,2]	= l1_data1[.,1]:^2
l1_data1[.,3]	= l1_data1[.,1]:^3
l1_data1[.,6]	= l1_data1[.,1]:*l1_data1[.,4]
l1_data1[.,7]	= l1_data1[.,1]:*l1_data1[.,5]
l1_data1[.,8]	= l1_data1[.,2]:*l1_data1[.,4]
l1_data1[.,9]	= l1_data1[.,2]:*l1_data1[.,5]
l1_data1[.,10]	= l1_data1[.,3]:*l1_data1[.,4]
l1_data1[.,11]	= l1_data1[.,3]:*l1_data1[.,5]
l1_data1[.,14]	= l1_data1[.,1]:*l1_data1[.,12]
l1_data1[.,15]	= l1_data1[.,1]:*l1_data1[.,13]
l1_data1[.,16]	= l1_data1[.,2]:*l1_data1[.,12]
l1_data1[.,17]	= l1_data1[.,2]:*l1_data1[.,13]
l1_data1[.,18]	= l1_data1[.,3]:*l1_data1[.,12]
l1_data1[.,19]	= l1_data1[.,3]:*l1_data1[.,13]

l1_data2=data2
l1_data2[.,1]	= l1_data2[.,1] - ones
l1_data2[.,2]	= l1_data2[.,1]:^2
l1_data2[.,3]	= l1_data2[.,1]:^3
l1_data2[.,6]	= l1_data2[.,1]:*l1_data2[.,4]
l1_data2[.,7]	= l1_data2[.,1]:*l1_data2[.,5]
l1_data2[.,8]	= l1_data2[.,2]:*l1_data2[.,4]
l1_data2[.,9]	= l1_data2[.,2]:*l1_data2[.,5]
l1_data2[.,10]	= l1_data2[.,3]:*l1_data2[.,4]
l1_data2[.,11]	= l1_data2[.,3]:*l1_data2[.,5]
l1_data2[.,14]	= l1_data2[.,1]:*l1_data2[.,12]
l1_data2[.,15]	= l1_data2[.,1]:*l1_data2[.,13]
l1_data2[.,16]	= l1_data2[.,2]:*l1_data2[.,12]
l1_data2[.,17]	= l1_data2[.,2]:*l1_data2[.,13]
l1_data2[.,18]	= l1_data2[.,3]:*l1_data2[.,12]
l1_data2[.,19]	= l1_data2[.,3]:*l1_data2[.,13]

while (i>0) {
	index	= data1[.,1]:<T+1				 
	i		= max(index)
	P		= normal(l1_data2*beta2p')
	s_hat 		= s_hat + exp((data3*beta3') - j*log(1+r)*ones):*(exp(l1_data1*beta1')):/(exp(l1_data1*beta1')+ exp(l1_data2*beta2'):*P):*index
	
	l1_data1=data1
	l1_data2=data2
	
	data1[.,1]	= data1[.,1] + ones
	data1[.,2]	= data1[.,1]:^2
	data1[.,3]	= data1[.,1]:^3
	data1[.,6]	= data1[.,1]:*data1[.,4]
	data1[.,7]	= data1[.,1]:*data1[.,5]
	data1[.,8]	= data1[.,2]:*data1[.,4]
	data1[.,9]	= data1[.,2]:*data1[.,5]
	data1[.,10]	= data1[.,3]:*data1[.,4]
	data1[.,11]	= data1[.,3]:*data1[.,5]
	data1[.,14]	= data1[.,1]:*data1[.,12]
	data1[.,15]	= data1[.,1]:*data1[.,13]
	data1[.,16]	= data1[.,2]:*data1[.,12]
	data1[.,17]	= data1[.,2]:*data1[.,13]
	data1[.,18]	= data1[.,3]:*data1[.,12]
	data1[.,19]	= data1[.,3]:*data1[.,13]

	data2[.,1]	= data2[.,1] + ones
	data2[.,2]	= data2[.,1]:^2
	data2[.,3]	= data2[.,1]:^3
	data2[.,6]	= data2[.,1]:*data2[.,4]
	data2[.,7]	= data2[.,1]:*data2[.,5]
	data2[.,8]	= data2[.,2]:*data2[.,4]
	data2[.,9]	= data2[.,2]:*data2[.,5]
	data2[.,10]	= data2[.,3]:*data2[.,4]
	data2[.,11]	= data2[.,3]:*data2[.,5]
	data2[.,14]	= data2[.,1]:*data2[.,12]
	data2[.,15]	= data2[.,1]:*data2[.,13]
	data2[.,16]	= data2[.,2]:*data2[.,12]
	data2[.,17]	= data2[.,2]:*data2[.,13]
	data2[.,18]	= data2[.,3]:*data2[.,12]
	data2[.,19]	= data2[.,3]:*data2[.,13]

	data3[.,1]	= data3[.,1] + ones
	data3[.,2]	= data3[.,1]:^2
	data3[.,3]	= data3[.,1]:^3
	data3[.,4]	= data3[.,4] + ones
	data3[.,5]	= data3[.,4]:^2
	data3[.,6]	= data3[.,4]:^3
	data3[.,9]	= data3[.,1]:*data3[.,7]
	data3[.,10]	= data3[.,1]:*data3[.,8]
	data3[.,11]	= data3[.,2]:*data3[.,7]
	data3[.,12]	= data3[.,2]:*data3[.,8]
	data3[.,13]	= data3[.,3]:*data3[.,7]
	data3[.,14]	= data3[.,3]:*data3[.,8]
	data3[.,17]	= data3[.,1]:*data3[.,15]
	data3[.,18]	= data3[.,1]:*data3[.,16]
	data3[.,19]	= data3[.,2]:*data3[.,15]
	data3[.,20]	= data3[.,2]:*data3[.,16]
	data3[.,21]	= data3[.,3]:*data3[.,15]
	data3[.,22]	= data3[.,3]:*data3[.,16]
	data3[.,70]	= data3[.,4]:*data3[.,68]
	data3[.,71]	= data3[.,4]:*data3[.,69]
	data3[.,72]	= data3[.,5]:*data3[.,68]
	data3[.,73]	= data3[.,5]:*data3[.,69]
	data3[.,74]	= data3[.,6]:*data3[.,68]
	data3[.,75]	= data3[.,6]:*data3[.,69]
	data3[.,78]	= data3[.,4]:*data3[.,76]
	data3[.,79]	= data3[.,4]:*data3[.,77]
	data3[.,80]	= data3[.,5]:*data3[.,76]
	data3[.,81]	= data3[.,5]:*data3[.,77]
	data3[.,82]	= data3[.,6]:*data3[.,76]
	data3[.,83]	= data3[.,6]:*data3[.,77]


	j = j+1
}

st_store(.,"s_hat",s_hat)
end

replace s_hat=s_hat/Q1
   
/* Define s and pi and q explicitly (s_hat already defined above) */ 

tsset person year
if $assets==1 {
	gen tot_assets3_temp=l2.tot_assets3
	drop tot_assets3
	ren tot_assets3_temp tot_assets3
}

egen med_assets = median(tot_assets3), by(age educ)
egen mean_assets = mean(tot_assets3), by(age educ) 

gen D3 	= tot_assets3
gen s 	= D1/(D1+D2)
gen pi 	= 1-(D3/(D1+D2+D3))
gen s_alt	= D1_alt/(D1_alt+D2_alt)
gen pi_alt	= 1-(D3/(D1_alt+D2_alt+D3))
gen pi_tax  = 1-(D3/(Q1+D3))

gen q=ly/(ly + wly)	/* Note that we will use the t-2 q in the implementation */ 


/*  ************************* NOTE *****************************
	There is a problem of missing values in consumption. One 
	approach would be to work from start only with non-missing
	values of consumption. Alternatively (what we do here) is 
	first estimate each process with as much observations as
	possible. This generates the strage effect that there are 
	more observations with missing consumption than observations
	with missing hours for wife. 
	It means that when keeping duyw!=. & duww!=. (see below)
	we still have missing observations for consumption. Note 
	however that assuming that the missing is random, this 
	should not be a problem. 
	***********************************************************/ 
	
*** In this version the residuals are already in growth rates
gen duc=uc
gen duy=uy
gen duyw=uyw
gen duw=uw
gen duww=uww
gen dutoty=utoty

gen l2_z_eta 	= l2.z_eta
gen l2_inv_mills= l2.inv_mills
gen l4_inv_mills= l4.inv_mills

*** For conditional means also need the difference of mills ratio
gen d_inv_mills = inv_mills-l2.inv_mills

/* generate unconditional moments */ 

keep person year duy duc duyw duw duww dutoty inv_mills z_eta l2_inv_mills l4_inv_mills d_inv_mills l2_z_eta s s_hat pi pi_tax q age agew yb fin_assets1 tot_assets3 house D1 D2 D3 log_c log_y log_yw log_h log_hw log_w log_ww med_assets mean_assets *educ* *_alt *help* tau

tsset person year

gen unc_duc2 		= duc^2
gen unc_duw2		= duw^2
gen unc_duy2 		= duy^2

gen unc_duy_lag		= duy*l2.duy
gen unc_duw_lag		= duw*l2.duw

gen unc_duw_duc		= duw*duc
gen unc_duw_duy		= duw*duy
gen unc_duc_duy		= duc*duy

gen unc_duc_lag		= duc*l2.duc /* Used only in the case of consumption measurement error */ 

gen unc_log_y		= log_y
gen unc_log_y2		= log_y^2
gen unc_log_yw		= log_yw
gen unc_log_yw2		= log_yw^2
gen unc_log_h		= log_h
gen unc_log_h2		= log_h^2

/*
xtile asset_pec = tot_assets3, n(4)
xi i.asset_pec, noomit
forvalues i=1/4 {
	ren _Iasset_pec_`i' asset_pec`i'
}
*/


/*xtile age_cat = age if duw!=. , n(5)
replace age_cat=4 if age_cat==5
*/
*xtile age_cat = age, n(5)
gen age_cat=1
replace age_cat=2 if age>=40
replace age_cat=3 if age>=45
replace age_cat=4 if age>=50
replace age_cat=5 if age>=55

/* 4 age categories */
*replace age_cat=4 if age_cat==5

/* 3 age categories */
*drop age_cat
*gen age_cat=1
*replace age_cat=2 if age>=40
*replace age_cat=3 if age>=55



tab age_cat, gen(age_cat)

/* Generate assets categories measures. Note that we always use the assets from the first year
   that the household is observed to generate the percentiles. */ 
egen min_year = min(year), by(person)
gen ass_temp = tot_assets3 
if $assets==1 {
	replace ass_temp = f2.tot_assets3
}
gen init_asset_temp = ass_temp if year==min_year
egen init_asset = min(init_asset_temp), by(person)
xtile asset_pec = init_asset, n(3)
drop init_asset* ass_temp

gen ass_temp = fin_assets1 
if $assets==1 {
	replace ass_temp = f2.fin_assets1 
}
gen init_asset_temp = ass_temp if year==min_year
egen init_asset = min(init_asset_temp), by(person)
xtile liqasset_pec = init_asset, n(3)
drop init_asset* ass_temp

xtile pi_pec = pi, n(3)

save estimation_input_unc`consumption'_rhs${rhs}_part${participation}, replace 

/* generate conditional moments */

keep if duyw!=. & duww!=. 			/* this captures our condition for P_t=1 and P_t-2=1 for the wife */ 
drop unc_*

/* Prepare all the cross residuals for the cross moments */ 

tsset person year

*** own and cross moments

gen duc2 		= duc^2
gen duw2		= duw^2
gen duww2		= duww^2
gen duy2 		= duy^2
gen duyw2 		= duyw^2

gen duy_lag		= duy*l2.duy
gen duyw_lag	= duyw*l2.duyw
gen duw_lag		= duw*l2.duw
gen duww_lag	= duww*l2.duww

gen duw_duc		= duw*duc
gen duw_duy		= duw*duy
gen duw_duyw	= duw*duyw
gen duww_duc	= duww*duc
gen duww_duy	= duww*duy
gen duww_duyw	= duww*duyw
gen duc_duy		= duc*duy
gen duc_duyw	= duc*duyw
gen duy_duyw	= duy*duyw

gen duw_duww	= duw*duww		/* 	note that one of our indentification assumptions is that this (and the following two) are zero 
									(zero will remain zero in the truncated case. Non-zero will need to be corrected */ 
gen duw_duww_lag= duw*l2.duww
gen duww_duw_lag= l2.duw*duww							
	
gen duc_lag		= duc*l2.duc 	/* Used only in the case of consumption measurement error */ 

gen duw_duc_lag = duw*l2.duc	/* The following 8 are only used for the nonseparable case of wages and hours */ 
gen duc_duw_lag = duc*l2.duw	
gen duy_duc_lag = duy*l2.duc
gen duc_duy_lag = duc*l2.duy
gen duww_duc_lag = duww*l2.duc	
gen duc_duww_lag = duc*l2.duww	
gen duyw_duc_lag = duyw*l2.duc
gen duc_duyw_lag = duc*l2.duyw

gen dutoty2 	= dutoty^2

*** Moments for Measurement error calculations

gen log_y2		= log_y^2		/* Used only for ME in income */
gen log_yw2		= log_yw^2
gen log_h2		= log_h^2

tab year, su(duc2)
tab year, su(duc_lag)

save estimation_input`consumption'_rhs${rhs}_part${participation}, replace 

/*  ************************* NOTE *****************************
	The number of observations going to the estimation for each
	of the moments is different. Other than the problem with 
	missing values for consumption, the lagged moments require
	have less observations by definition (since it only starts 
	2002 + it needs 4 the observation to be in 3 consecutive 
	waves and not only 2).
	***********************************************************/ 

*** Addition: Conditional moments will be defined by missing observations
u estimation_input_unc`consumption'_rhs${rhs}_part${participation}, clear
drop unc_*
gen constraint=(duyw!=. & duww!=.)

/* Prepare all the cross residuals for the cross moments */ 

tsset person year

*** own and cross moments
gen duc2 		= duc^2
gen duw2		= duw^2
gen duww2		= duww^2
gen duy2 		= duy^2
gen duyw2 		= duyw^2

gen duy_lag		= duy*l2.duy
gen duyw_lag	= duyw*l2.duyw
gen duw_lag		= duw*l2.duw
gen duww_lag	= duww*l2.duww

gen duw_duc		= duw*duc
gen duw_duy		= duw*duy
gen duw_duyw	= duw*duyw
gen duww_duc	= duww*duc
gen duww_duy	= duww*duy
gen duww_duyw	= duww*duyw
gen duc_duy		= duc*duy
gen duc_duyw	= duc*duyw
gen duy_duyw	= duy*duyw

gen duw_duww	= duw*duww		/* 	note that one of our indentification assumptions is that this (and the following two) are zero 
									(zero will remain zero in the truncated case. Non-zero will need to be corrected */ 
gen duw_duww_lag= duw*l2.duww
gen duww_duw_lag= l2.duw*duww							
	
gen duc_lag		= duc*l2.duc 	/* Used only in the case of consumption measurement error */ 

gen duw_duc_lag = duw*l2.duc	/* The following 16 are only used for the nonseparable case of wages and hours */ 
gen duc_duw_lag = duc*l2.duw	
gen duy_duc_lag = duy*l2.duc
gen duc_duy_lag = duc*l2.duy
gen duww_duc_lag = duww*l2.duc	
gen duc_duww_lag = duc*l2.duww	
gen duyw_duc_lag = duyw*l2.duc
gen duc_duyw_lag = duc*l2.duyw

gen duw_duy_lag = duw*l2.duy	
gen duy_duw_lag = duy*l2.duw	
gen duw_duyw_lag = duw*l2.duyw	
gen duyw_duw_lag = duyw*l2.duw	
gen duww_duy_lag = duww*l2.duy	
gen duy_duww_lag = duy*l2.duww	
gen duww_duyw_lag = duww*l2.duyw	
gen duyw_duww_lag = duyw*l2.duww	

gen dutoty2 	= dutoty^2

*** Moments for Measurement error calculations

gen log_y2		= log_y^2		/* Used only for ME in income */
gen log_yw2		= log_yw^2
gen log_h2		= log_h^2

tab year, su(duc2)
tab year, su(duc_lag)

sort person year
save estimation_input_all`consumption'_rhs${rhs}_part${participation}, replace 


