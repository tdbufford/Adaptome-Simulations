################################################################################
### RUN SIMULATIONS                                            #################
### Various Models                                             #################
################################################################################

library(R2jags)
library(dplyr)
library(tidyverse)
library(stringr)
library(coda)
library(MASS)
library(ebal)
library(ggplot2)


### RUN SIMULATIONS ##########################################################

source("PlatformTrial_Entropy.R")
source("PlatformTrial_EntropyCovars.R")
source("PlatformTrial_noEntropy.R")


# treatment effect for each adaptation
txEffects = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)

# variance for data generated from normal distribution
sigma2 = 1

# superiority threshold
supThresh = 0.9

# futility threshold
futThresh = 0.1

# number of treatment arms 
arms = 4 

# number of individuals per arm per interim analysis
NperInterim = 300

# maximum number of individuals total before stopping (all arms combined)
maxSample = 20000

### Sims with Entropy Balancing #############################################

# trial name
trialName = "EntropyBalancing_IncreasingEffects"

# run simulations for each outcome design
outcomeDesign = 1
for(i in 1:100){
  trialSimEntropy(txEffects, sigma2, arms, NperInterim, maxSample, 
               supThresh, futThresh, i, trialName, outcomeDesign)
}

outcomeDesign = 2
for(i in 1:100){
  trialSimEntropy(txEffects, sigma2, arms, NperInterim, maxSample, 
                  supThresh, futThresh, i, trialName, outcomeDesign)
}


# trial name
trialName = "EntropyBalancing_RandomEffects"

# run simulations for each outcome design
outcomeDesign = 1
for(i in 1:100){
  # shuffle order of tx effects
  txEffects_shuffle = sample(txEffects, size = length(txEffects), replace = FALSE)
  
  trialSimEntropy(txEffects_shuffle, sigma2, arms, NperInterim, maxSample, 
                  supThresh, futThresh, i, trialName, outcomeDesign)
}

outcomeDesign = 2
for(i in 1:100){
  # shuffle order of tx effects
  txEffects_shuffle = sample(txEffects, size = length(txEffects), replace = FALSE)
  
  trialSimEntropy(txEffects_shuffle, sigma2, arms, NperInterim, maxSample, 
                  supThresh, futThresh, i, trialName, outcomeDesign)
}

### Sims with Entropy Balancing + Covariates ###################################


# trial name
trialName = "EB_Covariates_IncreasingEffects"

# run simulations for each outcome design
outcomeDesign = 1
for(i in 1:100){
  trialSimEntropyC(txEffects, sigma2, arms, NperInterim, maxSample, 
                  supThresh, futThresh, i, trialName, outcomeDesign)
}

outcomeDesign = 2
for(i in 1:100){
  trialSimEntropyC(txEffects, sigma2, arms, NperInterim, maxSample, 
                  supThresh, futThresh, i, trialName, outcomeDesign)
}


# trial name
trialName = "EB_Covariates_RandomEffects"

# run simulations for each outcome design
outcomeDesign = 1
for(i in 1:100){
  # shuffle order of tx effects
  txEffects_shuffle = sample(txEffects, size = length(txEffects), replace = FALSE)
  
  trialSimEntropyC(txEffects_shuffle, sigma2, arms, NperInterim, maxSample, 
                  supThresh, futThresh, i, trialName, outcomeDesign)
}

outcomeDesign = 2
for(i in 1:100){
  # shuffle order of tx effects
  txEffects_shuffle = sample(txEffects, size = length(txEffects), replace = FALSE)
  
  trialSimEntropyC(txEffects_shuffle, sigma2, arms, NperInterim, maxSample, 
                  supThresh, futThresh, i, trialName, outcomeDesign)
}


### Sims with Neither (No EB, No Covariate Adjustment) #########################


# trial name
trialName = "Neither_IncreasingEffects"

# run simulations for each outcome design
outcomeDesign = 1
for(i in 1:100){
  trialSimNoEntropy(txEffects, sigma2, arms, NperInterim, maxSample, 
                   supThresh, futThresh, i, trialName, outcomeDesign)
}

outcomeDesign = 2
for(i in 1:100){
  trialSimNoEntropy(txEffects, sigma2, arms, NperInterim, maxSample, 
                   supThresh, futThresh, i, trialName, outcomeDesign)
}


# trial name
trialName = "Neither_RandomEffects"

# run simulations for each outcome design
outcomeDesign = 1
for(i in 1:100){
  # shuffle order of tx effects
  txEffects_shuffle = sample(txEffects, size = length(txEffects), replace = FALSE)
  
  trialSimNoEntropy(txEffects_shuffle, sigma2, arms, NperInterim, maxSample, 
                  supThresh, futThresh, i, trialName, outcomeDesign)
}

outcomeDesign = 2
for(i in 1:100){
  # shuffle order of tx effects
  txEffects_shuffle = sample(txEffects, size = length(txEffects), replace = FALSE)
  
  trialSimNoEntropy(txEffects_shuffle, sigma2, arms, NperInterim, maxSample, 
                  supThresh, futThresh, i, trialName, outcomeDesign)
}

