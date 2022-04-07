################################################################################
### RUN SIMULATIONS                                            #################
### BASED ON SCHOOL DATA                                       #################
################################################################################


txEffects = c(0, -1.9, 0.07, 1.4, 0.4, 0.85, 0.72, 0.11, 0.77)
sigma2 = 1
arms = 4
NperInterim = 100
maxSample = 2000
supThresh = 0.8
futThresh = 0.2
trialNum = 1
trialName = "SchoolData"
source("PlatformTrialSim_SchoolData.R")

NperInterim = 40
for(i in 1:100){
  trialSimEntropySchool(txEffects, sigma2, arms, NperInterim, maxSample,
                        supThresh, futThresh, i, trialName)
}

NperInterim = 60
for(i in 1:100){
  trialSimEntropySchool(txEffects, sigma2, arms, NperInterim, maxSample,
                        supThresh, futThresh, i, trialName)
}

NperInterim = 80
for(i in 1:100){
  trialSimEntropySchool(txEffects, sigma2, arms, NperInterim, maxSample,
                        supThresh, futThresh, i, trialName)
}

NperInterim = 100
for(i in 1:100){
  trialSimEntropySchool(txEffects, sigma2, arms, NperInterim, maxSample,
                        supThresh, futThresh, i, trialName)
}

NperInterim = 120
for(i in 1:100){
  trialSimEntropySchool(txEffects, sigma2, arms, NperInterim, maxSample,
                        supThresh, futThresh, i, trialName)
}

NperInterim = 140
for(i in 1:100){
  trialSimEntropySchool(txEffects, sigma2, arms, NperInterim, maxSample,
                        supThresh, futThresh, i, trialName)
}

NperInterim = 160
for(i in 1:100){
  trialSimEntropySchool(txEffects, sigma2, arms, NperInterim, maxSample,
                        supThresh, futThresh, i, trialName)
}

NperInterim = 180
for(i in 1:100){
  trialSimEntropySchool(txEffects, sigma2, arms, NperInterim, maxSample,
                        supThresh, futThresh, i, trialName)
}

NperInterim = 200
for(i in 1:100){
  trialSimEntropySchool(txEffects, sigma2, arms, NperInterim, maxSample,
                        supThresh, futThresh, i, trialName)
}
