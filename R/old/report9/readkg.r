#log experience profile
kappa = read.table("old/Kaplan\ and\ Giovanni/Input/kappasmooth.txt")$V1 #25~59
kappa = kappa + log(0.75)

#survival probabilities
surprob = read.table("old/Kaplan\ and\ Giovanni/Input/surprobsmooth.txt")$V1 #60~89

#age-specific variances
varetaage = read.table("old/Kaplan\ and\ Giovanni/Input/varetaage.txt")$V1 #26~59
Vetavec = 0.01*varetaage

#unconsurprob annprem
unconsurprob = rep(1,35)
for(it in 2:35){
  unconsurprob[it] = unconsurprob[it-1]* surprob[it-1]
}

initwealthdist=read.table("old/Kaplan\ and\ Giovanni/Input/initwealthdist.txt")