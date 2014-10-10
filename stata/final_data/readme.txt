"data.dta" is the pre-selection data
"first_stage_regs.dta" is the result sample. muc, mutoty​, mua are the residues for consumption, income, and asset without moving deterministic trend. 

To get the result:
1. Run "select.do", which reads "data.dta", and generates "select.dta".
2. Run "first_stage_regs.do", which reads "select.dta", and generates "first_stage_regs.dta".
