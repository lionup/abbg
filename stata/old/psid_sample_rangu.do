******DOES SAMPLE SELECTION *************

/* 	NOTE: The concept has changed from the AER sampling:
	1. 	In this version I drop "vet" since need to update that from years before 99. 
	2. 	Max(educ) for wife is problematic unless checking that wife did not change. 
	3.	For this estimation procedure we need only pairs of years, therefore there are no restrictions on the number of consecutive years of non-missing data. 
		If there is any break (family change and so on) the person is brough back with new id. 
	4. 	Note that family change (fchg) is missing only for 93-98 in our sample. Since our sample starts at 98, and changes are not checked for the first year, 
		this is ok. 
	5. For now I keep the income sampling by the growth in household total income and not by growth in wages/earnings for head/wife
	Technically - I do it by dropping individual observations (and not the whole household) based on the 
	sampling criteria (e.g. not married), and in the end, whenever there is a household with missing 
	years in the data, I assing the years after the first missing with a new id. 
	6. Note that  topcoding is not an issue (only 2 top coded observations in the cleaned sample) for both constumption income and assets) 
	*/
	
*cd "$output"
clear all
set mem 100m
set more off
cap log close
log using notsampled, replace 

u psid_cons_notsampled,clear
drop pid

/* update info which was not asked in all years */ 
tsset id year
* replace vet = L1.vet if newhd==5 & year>=94 & vet==.
drop vet
gen avhy=ly/hours 					/* 	avhy is not reported by PSID starting the files of 1993. 
									    Before it is topcoded at 99.99 therefore I recalculte for all years */ 

* Makes the variable race consistent over the years
* Now 1 is white, 2 black, 3 others
forvalues i=1/4 {
	gen     rc=race`i' if race`i'<=2
	replace rc=3    if race`i'>=3 & race`i'<=7
	drop race`i'
	ren rc race`i'

	gen     wrc=wrace`i' if wrace`i'<=2
	replace wrc=3    if wrace`i'>=3 & wrace`i'<=7
	drop wrace`i'
	ren wrc wrace`i'
}

/*
Note: Reported switch in race across panel years is common. At least one reason is that there 
----  is an option to mention more than one race, and the respondants reverse the order of the mentions.
      This is especially true when there is a slight change in the categories. In 2005 survey for example, 
	  the 3rd category was change from "Native American" to "American Indian or Alaska Native". We deal 
	  with that by assuming that the clearest (and most stable) categories are white (1) and african-
	  american(2). Whenever a respondant mentions in any year an additional category he is assigned to the 
	  third category (other). 

*/

gen race=race1
replace race=3 if race2!=.		/* Any case of more than one reported race is also considered "other" */ 
egen race_temp = max(race), by(person)
replace race=race_temp

gen wrace=wrace1
replace wrace=3 if wrace2!=.
egen wrace_temp = max(wrace), by(person)
replace wrace = wrace_temp

drop *temp wrace1 wrace2 wrace3 wrace4 race1 race2 race3 race4

/*
Replaces educ with educ from individual file when missing. 
Makes the variable education consistent over the years
Now 1 is 0-11 grades (includes those with no schooling, i.e. can't read/write people); 
    2 is 12 grades, i.e. high school w or w/o nonacademic training;
    3 is college dropout, college degree, or collage & advanced/professional degree;
    a missing point denotes na/dk.
*/

replace educ_ind =. if educ_ind==99
replace weduc_ind=. if weduc_ind==99
replace educ = educ_ind if educ==.
replace weduc = weduc_ind if weduc==.

gen sc=1 if educ>=0  & educ<=11           	/* 0-11  grades*/
replace sc=2 if educ==12                        /* High school or 12 grades+nonacademic training*/
replace sc=3 if educ>=13 & educ<=17          	/* College dropout, BA degree, or collage & adv./prof. degree*/
replace sc=. if educ>17                    		/* Missing, NA, DK*/
drop educ
ren sc educ

gen sc=1 if weduc>=0  & weduc<=11           /* 0-11  grades*/
replace sc=2 if weduc==12                       /* High school or 12 grades+nonacademic training*/
replace sc=3 if weduc>=13 & weduc<=17          	/* College dropout, BA degree, or collage & adv./prof. degree*/
replace sc=. if weduc>17                    	/* Missing, NA, DK*/
drop weduc
ren sc weduc

egen maxed=max(educ),by(person)
gen educ2=educ
replace educ=maxed
drop maxed educ2

egen maxed=max(weduc),by(person)
gen educ2=weduc
replace weduc=maxed
drop maxed educ2 


* Design a demographically unstable household as one where some family composition change (apart from changes in
* people other than head-wife takes place).
* Then drop the years after the change for households with unstable demographical pattern

/*
fchg=0 (No change), 1 (Change in members other than Head or Wife), 2 (Head same, wife left,die, or is new)
3 (wife is now head), 4 (ex female head married and huisband is now head), 5 (some sample member other than 
ex Head or Wife is now head)
7 (ex wife head because husband in inst., now husband back and is head), 8 (Other) 
*/

******PAY ATTENTION TO MARITAL STATUS ETC
sort person year
egen miny 		= min(year),by(person)
gen break_d 	= (fchg>1 & year>98)
replace break_d	= 0 if year==miny & break_d==1
drop if break_d == 1
drop miny break_d

/* dropped the non-married heads */ 
drop if marit!=1 | agew==0

/* drop if state is missing */ 
drop if state==.|state==0|state==99

# delimit;
gen     region=1 if state==6  | state==18 | state==20 | state==28 | state==29 | state==31 
                              | state==37 | state==38 | state==44; 						 /* North East*/
replace region=2 if state==12 | state==13 | state==14 | state==15 | state==21 | state==22
                              | state==24 | state==26 | state==33 | state==34 | state==40 | state==48; /*Midwest*/
replace region=3 if state==1  | state==3  | state==7  | state==8  | state==9  | state==10 | state==16  
                              | state==17 | state==19 | state==23 | state==32 | state==35 
                              | state==39 | state==41 | state==42 | state==45 | state==47; 		 /*South*/
replace region=4 if state==2  | state==4  | state==5  | state==11 | state==25 | state==27 | state==30 
                              | state==36 | state==43 | state==46 | state==49 | state==50 | state==51; /*West*/                              
#delimit cr

/* Drop female household heads*/
drop if sex==2

/* Recode age so that there is no gap or jump*/
/* FIrst drop those with always missing info on age */

egen n=sum(person!=.),by(person)
egen na=sum(age==.),by(person)
drop if n==na
drop n na

egen n=sum(person!=.),by(person)
egen na=sum(agew==.),by(person)
drop if n==na
drop n na

egen lasty=max(year) if age!=.,by(person)
gen lastage=age if year==lasty
gen b=year-lastage
replace b=0 if b==.
egen yb=sum(b),by(person)
replace age=year-yb
drop lasty lastage b

egen lasty=max(year) if agew!=.,by(person)
gen lastage=agew if year==lasty
gen b=year-lastage
replace b=0 if b==.
egen wyb=sum(b),by(person)
replace agew=year-wyb
drop lasty lastage b

* Takes into account the retrospective nature of the data
replace age=age-1
replace agew=agew-1

/*

/* Drop additional extreme outlyers - identified from J Laird Summer 09*/

drop if out==1 | out==2 						/* drop additional outliers */ 

*/

/*  For intermittent "headship" drop the years after the break for long_sample=0 and bring back the family with new id one year after the break for long_sample=1 */
/* 	Also dropping if less than 2 consecutive periods in the sample */ 
/* 	This was moved to the end since intermittent "headship" might be due to missing race / educ etc. This way we are not dropping the whole
	HH but just the years after the missing data. Alternatively we can not drop these at all, but change the mindist_AER to account for "holes"
	in the series */ 

sort person year
qby person: gen dyear	= year-year[_n-1]
qby person: gen break_d = (dyear>1 & dyear!=. & year<=96 | dyear>2 & dyear!=. & year>96)
qby person: gen b_year 	= year if break_d == 1
egen break_year 		= min(b_year), by(person)

gen long_sample=0 			/* long_sample is equal to one for families which were broken and than put back in to the sample */ 
local ind=1
while `ind'>0 {
	sum person
	local max_id	=	r(max)
	count if year==break_year   /* for debugging */
	di r(N)						/* for debugging */ 
	* drop if year==break_year
	replace person = person + `max_id' if year>=break_year
	replace long_sample = 1 if year>break_year
	drop dyear break_d b_year break_year
	sort person year
	qby person: gen dyear	= year-year[_n-1]
	qby person: gen break_d = (dyear>1 & dyear!=. & year<=96 | dyear>2 & dyear!=. & year>96)
	qby person: gen b_year 	= year if break_d == 1
	egen break_year 		= min(b_year), by(person)
	count if break_d == 1
	local ind=r(N)
	di `ind'					/* for debugging */
	}


	/* merge minimum wage data */ 
sort state
merge state using psid_states
drop _merge

sort state_st year
merge state_st year using min_wage
drop if _merge!=3
drop _merge

*** merge prices
sort year
merge m:1 year using price_indices
keep if _m==3
drop _m

/* counting the number of observation for baseline sample */
drop if age<25 | age>65 
 
log close
saveold data3_2, replace 




