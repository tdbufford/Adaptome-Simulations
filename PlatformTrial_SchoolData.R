###############################################################################
## Platform Clinical Trial Simulation                                     ####
### Normal outcomes                                                       ####
### with Entropy Balancing                                                ####
### BASED ON SCHOOL DATA
###############################################################################

library(tidyverse)
library(dplyr)
library(R2jags)
library(stringr)
library(coda)
library(MASS)
library(ebal)
library(ggplot2)

source("helper_checkModelConvergenceFunction.R")
source("helper_modelDiagnosticGraphs.R")
source("helper_trialSummaryFunction.R")
source("helper_trialVarianceSummaryFunction.R")

## trialSimEntropySchool function #################################################
# 
# Simulates a Bayesian platform clinical trial.
# Treatment outcomes are continuous (assumed normally distributed).
# Treatments are dropped for futility and new ones are added along the way.
# At each interim analysis a model is fit to estimate mean outcomes for
# each treatment group. Higher is assumed to be better.
# The model is fit using JAGS to get posterior estimates
# 
# @param txEffects - vector of true means for each treatment
# @param sigma2 - the variance of the simulated data
# @param arms - number of active treatment arms - NOTE: code does not currently 
#                   accomodate more than 4. To increase, edit code at line 459
# @param NperInterim - number of participants (per group) added for each interim 
#                     analysis
# @param maxSample - maximum total sample size allowed (trial will stop when 
#                     limit reached)
# @param supThresh - superiority threshold, the probability of a treatment 
#                     being superior to the current reference group must meet 
#                     this threshold in order for the new treatment to be 
#                     switched to become the new reference group
# @param futThresh - when the probability that a treatment is better than the 
#                     current reference group falls below the futility threshold
#                     then the treatment arm is dropped and a new one is added
# @param trialNum - counter to keep track when doing multiple simulations
# @param trialName - name of trial, also for keeping track with multiple simulation 
#              versions
#              
# @return: no object is returned by the function, however several files are
#          created and saved. 
#          1) platformTrialSummary_'trialName'.csv - a summary of estimates and 
#             superiority probabilities for each interim analysis
#          2) totalSampleSizes_'trialName'.csv - a list of the total sample size 
#             for each trial (useful if running multiple simulations)
#          3) ModelConvergence_'trialName'.csv - a log of any model convergence 
#             errors that occurred (model is re-fit with longer chains and more
#             thinning if this occurs)
# 

trialSimEntropySchool = function(txEffects, sigma2, arms, NperInterim, maxSample, 
                        supThresh, futThresh, trialNum, trialName){
  
  
  ## 1) Initialize Variables
  
  # counter for number of adaptations that have been tested
  # number of arms is number of adaptations being tested +1 for control
  numAdapts = arms
  Arms = arms
  
  # For each platform trial, shuffle the labels/order of the tx groups
  # that way covariate imbalances are not always appearing in the same order
  order_grps = seq(1:length(txEffects))
  order_grps = sample(order_grps, size = length(txEffects), replace = FALSE)
  
  # Re-order treatment effects (same order as tx assignment based on covariates)
  # Note: DO NOT shuffle order of txEffects BEFORE inputting to the function
  txEffects = cbind(txEffects[order_grps[1]],
                    txEffects[order_grps[2]],
                    txEffects[order_grps[3]],
                    txEffects[order_grps[4]],
                    txEffects[order_grps[5]],
                    txEffects[order_grps[6]],
                    txEffects[order_grps[7]],
                    txEffects[order_grps[8]],
                    txEffects[order_grps[9]])
  
  # true parameters
  theta = txEffects[1:arms]
  
  
  # names of treatments being tested
  txNames = c("Adaptation 1 (Ref)")
  for(i in 2:numAdapts){
    txNames = c(txNames, paste("Adaptation", i, sep = " "))
  }
  
  # corresponding treatment numbers
  txNums = seq(from = 1, to = numAdapts, by = 1)
  
  # flag for which group is reference group
  refFlag = c(1, rep(0, numAdapts - 1))
  
  # initialize counter for number of interim analyses
  k = 1
  
  # initialize total sample size counter
  totalSample = arms * NperInterim
  
  # initialize vector to keep track of sample size per adaptation
  sampleSize = rep(0, arms)
  
  # initialize model convergence flag
  modelConverge = TRUE
  
  # counter for whether maximum number of iterations has been reached
  # want to do max of 3 interim analyses after last adaptation added
  maxIteration = -1
  
  
  
  ## 2) Initialize Priors for Bayesian Analysis
  # the priors are conjugate priors

  # 1 row per treatment arm
  # col1 is mean, and col2 is var for Normal(mean, var) prior
  # initialize with vague prior Normal(0, 10)
  thetapriors = cbind(rep(0, arms), rep(100, arms))
  
  # variance ~ Inverse Gamma(a', b')
  # precesion = 1/variance
  # precision ~ Gamma(a, b)
  sigmaprior = c(2, 2)
  
  
  ## Create external sample to balance against for entropy balancing
  
  #Simulate covariates based on LAUSD data
  #athletic eligibility GPA - use a truncated t-distribution
  x1 = rt(40000, df = 1, ncp = 0) + 2.5
  x1 = x1[which(x1 <= 4.0 & x1 >= 0)]
  x1 = x1[1:10000]
  
  #gender
  x2 = rbinom(10000, size = 1, prob = 0.53)
  
  #ptsd screen
  x3 = rbinom(10000, size = 1, prob = 0.22)
  
  #special education
  x4 = rbinom(10000, size = 1, prob = 0.09)
  
  external_sample = cbind(x1, x2, x3, x4)
  external_sample = getsquares(external_sample)
  
  
  #Betas estimated from logit models fit on real school data
  beta2 = as.matrix(c(-0.123, 0.048, 1.294, -0.278, 1, 0, 0, 0, 0, 0, 0, 0))
  beta3 = as.matrix(c(-0.054, -0.359, 1.115, 0.578, 0, 1, 0, 0, 0, 0, 0, 0))
  beta4 = as.matrix(c(-0.290, -0.251, 0.544, -0.342, 0, 0, 1, 0, 0, 0, 0, 0))
  beta5 = as.matrix(c(-0.811, 0.026, 0.917, -0.107, 0, 0, 0, 1, 0, 0, 0, 0))
  beta6 = as.matrix(c(-0.040, -0.336, 0.476, -0.299, 0, 0, 0, 0, 1, 0, 0, 0))
  beta7 = as.matrix(c(-0.364, -0.253, 1.106, -0.824, 0, 0, 0, 0, 0, 1, 0, 0))
  beta8 = as.matrix(c(-0.563, -0.253, 1.131, 0.262, 0, 0, 0, 0, 0, 0, 1, 0))
  beta9 = as.matrix(c(0.159, -0.202, 1.152, -1.118, 0, 0, 0, 0, 0, 0, 0, 1))
  
  
  ## 3) Loop for Platform Trial
  
  # continue simulation until max sample size is reached or 
  # model fails to converge 
  while(totalSample <= maxSample && modelConverge == TRUE){

    ## 3.1) Update Counters
    
    # number of newly recruited participants
    N = arms * NperInterim
    
    # update sample size tracker (sample size per adaptation)
    sampleSize = sampleSize + rep(NperInterim, arms)
    
    
    ## 3.2) Generate Data
    datasize = length(txEffects) * NperInterim * 3
    print(paste("length txEffects = ", length(txEffects)))
    print(paste("Nperinterim = ", NperInterim))
    print(paste("datasize = ", datasize))
    
    #athletic eligibility GPA - use a truncated t-distribution
    x1 = rt(datasize*4, df = 1, ncp = 0) + 2.5
    x1 = x1[which(x1 <= 4.0 & x1 >= 0)]
    x1 = x1[1:datasize]
    #gender
    x2 = rbinom(datasize, size = 1, prob = 0.53)
    #ptsd screen
    x3 = rbinom(datasize, size = 1, prob = 0.22)
    #special education
    x4 = rbinom(datasize, size = 1, prob = 0.09)
    # x4 = rbinom(datasize, size = 1, prob = 0.4)
    
    epsilon1 = rnorm(datasize, mean = 0, sd = 4)
    epsilon2 = rnorm(datasize, mean = 0, sd = 4)
    epsilon3 = rnorm(datasize, mean = 0, sd = 4)
    epsilon4 = rnorm(datasize, mean = 0, sd = 4)
    epsilon5 = rnorm(datasize, mean = 0, sd = 4)
    epsilon6 = rnorm(datasize, mean = 0, sd = 4)
    epsilon7 = rnorm(datasize, mean = 0, sd = 4)
    epsilon8 = rnorm(datasize, mean = 0, sd = 4)

    X = cbind(x1, x2, x3, x4, epsilon1, epsilon2, epsilon3, epsilon4, 
            epsilon5, epsilon6, epsilon7, epsilon8)
    
    
    # Multinomial logit model for determining treatment groups
    # using a set of independent binary regressions 
    # R is the number of possible treatments (ie. length of txEffects)
    
    # ln[Pr(T = 2)/Pr(T = 1)] = Beta2.1*X[, 1] + Beta2.2*X[, 2] + Beta2.3*X[, 3] + Beta2.4*X[, 4] + Beta2.5*X[, 5] + Beta2.6*X[, 6]
    # .
    # .
    # ln[Pr(T = R)/Pr(T = 1)] = BetaR.1*X[, 1] + BetaR.2*X[, 2] + BetaR.3*X[, 3] + BetaR.4*X[, 4] + BetaR.5*X[, 5] + BetaR.6*X[, 6]
    
    # Probabilities of y_i ending up in each group:

    #Pr(T = 1) = 1 / [1 + exp(Beta2*x_i) + ... + exp(BetaR * x_i)]
    #Pr(T = 2) = exp(Beta2*x_i) / [1 + exp(Beta2*x_i) + ... + exp(BetaR * x_i)]
    # .
    # .
    #Pr(T = R) = exp(BetaR*x_i) / [1 + exp(Beta2*x_i) + ... + exp(BetaR * x_i)]
    
    
    # LINEAR MODEL FOR ASSIGNMENT OF OBSERVATIONS TO GROUPS
    score2 = X %*% beta2
    score3 = X %*% beta3
    score4 = X %*% beta4
    score5 = X %*% beta5
    score6 = X %*% beta6
    score7 = X %*% beta7
    score8 = X %*% beta8
    score9 = X %*% beta9
    
    denominator = 5 + exp(score2) + exp(score3) + exp(score4) + exp(score5) + 
      exp(score6) + exp(score7) + exp(score8) + exp(score9) 
    
    prT1 = 5/denominator
    prT2 = exp(score2)/denominator
    prT3 = exp(score3)/denominator
    prT4 = exp(score4)/denominator
    prT5 = exp(score5)/denominator
    prT6 = exp(score6)/denominator
    prT7 = exp(score7)/denominator
    prT8 = exp(score8)/denominator
    prT9 = exp(score9)/denominator
    
    probs = cbind(prT1, prT2, prT3, prT4, prT5, prT6, prT7, prT8, prT9)
    
    
    # calculate probabilities of ending up in each tx group (for each observation)
    # based on multinomial logit model
    probs = cbind(probs[, order_grps[1]],
                  probs[, order_grps[2]],
                  probs[, order_grps[3]],
                  probs[, order_grps[4]],
                  probs[, order_grps[5]],
                  probs[, order_grps[6]],
                  probs[, order_grps[7]],
                  probs[, order_grps[8]],
                  probs[, order_grps[9]])
    
    # assign treatment groups based on probabilities, which are based on covariates
    Treatment = apply(probs, 1, function(x) rmultinom(1, 1, x))
    Treatment = as.data.frame(t(Treatment)) %>%
      rename(t1 = V1, 
             t2 = V2,
             t3 = V3,
             t4 = V4,
             t5 = V5,
             t6 = V6,
             t7 = V7,
             t8 = V8,
             t9 = V9) %>%
      mutate(treatment = case_when(
        t1 == 1 ~ 1,
        t2 == 1 ~ 2,
        t3 == 1 ~ 3,
        t4 == 1 ~ 4,
        t5 == 1 ~ 5,
        t6 == 1 ~ 6,
        t7 == 1 ~ 7,
        t8 == 1 ~ 8,
        t9 == 1 ~ 9)) -> Treatment
    
    # generate random variation
    eta = rnorm(datasize, mean = 0.3, sd = 3.7)
    
    #Relationship between outcome and covariates estimated from real LAUSD data
    # #Outcome Design: 
    Outcome = 
      data.frame(
        outcome = (0.57*X[, 1] + 0.3*X[, 2] - 0.89*X[, 3] + 0.08*X[, 4] + 8.7 + eta))

    # combine treatment indicators and covariates into data frame
    Treatment %>%
      mutate(x1 = X[, 1],
             x2 = X[, 2],
             x3 = X[, 3],
             x4 = X[, 4]) -> data
    
    # add outcome to data frame
    data <- as.data.frame(cbind(data, Outcome))
    
    
    # check to make sure enough participants in each group
    data %>% 
      group_by(treatment) %>% 
      summarise(n = n()) %>% 
      mutate(n_flag = ifelse(n < NperInterim, 1, 0)) %>%
      summarise(need_more_data = sum(n_flag)) -> need_more_data
    need_more_data = need_more_data[[1]]
    
    # if not enough, generate more data
    if(need_more_data > 0){

      print("need more data")

      #athletic eligibility GPA - use a truncated t-distribution
      x1 = rt(datasize*4, df = 1, ncp = 0) + 2.5
      x1 = x1[which(x1 <= 4.0 & x1 >= 0)]
      x1 = x1[1:datasize]
      #gender
      x2 = rbinom(datasize, size = 1, prob = 0.53)
      #ptsd screen
      x3 = rbinom(datasize, size = 1, prob = 0.22)
      #special education
      x4 = rbinom(datasize, size = 1, prob = 0.09)
      # x4 = rbinom(datasize, size = 1, prob = 0.4)

      epsilon1 = rnorm(datasize, mean = 0, sd = 4)
      epsilon2 = rnorm(datasize, mean = 0, sd = 4)
      epsilon3 = rnorm(datasize, mean = 0, sd = 4)
      epsilon4 = rnorm(datasize, mean = 0, sd = 4)
      epsilon5 = rnorm(datasize, mean = 0, sd = 4)
      epsilon6 = rnorm(datasize, mean = 0, sd = 4)
      epsilon7 = rnorm(datasize, mean = 0, sd = 4)
      epsilon8 = rnorm(datasize, mean = 0, sd = 4)

      X = cbind(x1, x2, x3, x4, epsilon1, epsilon2, epsilon3, epsilon4,
                epsilon5, epsilon6, epsilon7, epsilon8)

      score2 = X %*% beta2
      score3 = X %*% beta3
      score4 = X %*% beta4
      score5 = X %*% beta5
      score6 = X %*% beta6
      score7 = X %*% beta7
      score8 = X %*% beta8
      score9 = X %*% beta9

      denominator = 5 + exp(score2) + exp(score3) + exp(score4) + exp(score5) +
        exp(score6) + exp(score7) + exp(score8) + exp(score9)

      prT1 = 5/denominator
      prT2 = exp(score2)/denominator
      prT3 = exp(score3)/denominator
      prT4 = exp(score4)/denominator
      prT5 = exp(score5)/denominator
      prT6 = exp(score6)/denominator
      prT7 = exp(score7)/denominator
      prT8 = exp(score8)/denominator
      prT9 = exp(score9)/denominator

      probs = cbind(prT1, prT2, prT3, prT4, prT5, prT6, prT7, prT8, prT9)

      probs = cbind(probs[, order_grps[1]],
                    probs[, order_grps[2]],
                    probs[, order_grps[3]],
                    probs[, order_grps[4]],
                    probs[, order_grps[5]],
                    probs[, order_grps[6]],
                    probs[, order_grps[7]],
                    probs[, order_grps[8]],
                    probs[, order_grps[9]])

      Treatment = apply(probs, 1, function(x) rmultinom(1, 1, x))
      Treatment = as.data.frame(t(Treatment)) %>%
        rename(t1 = V1,
               t2 = V2,
               t3 = V3,
               t4 = V4,
               t5 = V5,
               t6 = V6,
               t7 = V7,
               t8 = V8,
               t9 = V9) %>%
        mutate(treatment = case_when(
          t1 == 1 ~ 1,
          t2 == 1 ~ 2,
          t3 == 1 ~ 3,
          t4 == 1 ~ 4,
          t5 == 1 ~ 5,
          t6 == 1 ~ 6,
          t7 == 1 ~ 7,
          t8 == 1 ~ 8,
          t9 == 1 ~ 9)) -> Treatment

      eta = rnorm(datasize, mean = 0.3, sd = 3.7)
      
      #Relationship between outcome and covariates estimated from real LAUSD data
      # #Outcome Design: 
      Outcome = 
        data.frame(
          outcome = (0.57*X[, 1] + 0.3*X[, 2] - 0.89*(X[, 3]>=3) + 0.08*X[, 4] + 8.7 + eta))

      # combine treatment indicators and covariates into data frame
      Treatment %>%
        mutate(x1 = X[, 1],
               x2 = X[, 2],
               x3 = X[, 3],
               x4 = X[, 4]) -> data_supplement

      data_supplement <- as.data.frame(cbind(data_supplement, Outcome))

      data = as.data.frame(rbind(data, data_supplement))
    }
    
    
    # Add in treatment effects
    for(i in 1:length(txEffects)){
      data %>%
        mutate(outcome = outcome +
                 txEffects[i]*(treatment == i)) -> data
    }
    
    # restrict range of outcome from 0 to 15
    data %>%
      mutate(outcome = ifelse(outcome < 0, 0, outcome), 
             outcome = ifelse(outcome > 15, 15, outcome)) -> data
    
    # increase number of observations with outcome 0, 5, 10, or 15
    random = rbinom((datasize + need_more_data*datasize), 1, prob = 0.20)
    data %>%
      mutate(chance = random) %>%
      mutate(outcome = ifelse(outcome < 5 & chance == 1, 0, outcome), 
             outcome = ifelse(outcome > 5 & outcome < 7.5 & chance == 1, 5, outcome),
             outcome = ifelse(outcome >= 7.5 & outcome < 15 & chance == 1, 10, outcome)
             ) -> data
    
    # check distribution of outcome 
    hist(data$outcome)
    mean(data$outcome)
    median(data$outcome)
    var(data$outcome)
    
    # restrict dataset to be analyzed to the active treatment arms 
    # and specified number of observations per treatment group
    # !!! Need to add another if statement if want more than 4 active treatment arms
    if(arms == 4){
      data %>% filter(treatment == txNums[1] |
                      treatment == txNums[2] |
                      treatment == txNums[3] |
                      treatment == txNums[4]) %>%
               group_by(treatment) %>%
               slice_head(n = NperInterim) %>%
               ungroup() -> dataset
      
    }else if(arms == 3){
      data %>% filter(treatment == txNums[1] |
                        treatment == txNums[2] |
                        treatment == txNums[3]) %>%
        group_by(treatment) %>%
        slice_head(n = NperInterim) %>%
        ungroup() -> dataset
      
    }else if(arms == 2){
      data %>% filter(treatment == txNums[1] |
                      treatment == txNums[2]) %>%
        group_by(treatment) %>%
        slice_head(n = NperInterim) %>%
        ungroup() -> dataset
      
    }else if(arms == 1){
      data %>% filter(treatment == txNums[1]) %>%
        group_by(treatment) %>%
        slice_head(n = NperInterim) %>%
        ungroup() -> dataset
    }
    
    
    ### Check covariate overlap #############################################
    # Optional code: uncomment this if you want to see some graphs of the generated data 
    # will need to add break points to pause at graphs
    #
    # data %>%
    #   ggplot() +
    #   geom_density(aes(x = x1, group = treatment, color = as.factor(treatment))) +
    #   theme_bw()
    # 
    # data %>%
    #   ggplot() +
    #   geom_density(aes(x = x2, group = treatment, color = as.factor(treatment))) +
    #   theme_bw()
    # 
    # data %>%
    #   ggplot() +
    #   geom_density(aes(x = x3, group = treatment, color = as.factor(treatment))) +
    #   theme_bw()
    # 
    # data %>%
    #   ggplot() +
    #   geom_density(aes(x = x4, group = treatment, color = as.factor(treatment))) +
    #   theme_bw()
    #############################################################################
    
    ### APPLY ENTROPY BALANCING FOR ACTIVE TREATMENTS
    
    # To estimate ATE (Average Treatment Effect) we want to assess covariate balance 
    # between a given treatment group and the entire population
    
    # By default, the entropy balancing function assumes 0 is the control group 
    # and 1 is the treatment group. The treatment group is used to calculate 
    # target moments for covariate distributions, and control group is 
    # weighted to match
    
    # We want to weight each treatment group to match covariate moments of entire
    # population, so we need entire population with 1 as treatment indicator plus
    # repeated observations from treatment group of interest with 0 as indicator
    
    # We perform entropy balancing separately for each treatment group to weight it
    # to match the entire population
    
    
    # square covariates to use as input to entropy balancing
    # this allows us to match the variances for covariate distributions
    X_squares = getsquares(dataset[,11:14])
    X_squares = data.frame(X_squares)
    
    cbind(X_squares, dataset[,-(11:14)]) -> dataset
    
    
    # select each tx group and set treatment indicator to 0
    # save each as a separate data frame called "txgroup*"
    for(i in 1:arms){
      assign(paste("txgroup", i, sep=""), 
             (dataset %>% 
                     filter(treatment == txNums[i]) %>% 
                     mutate(treatment = 0)
              )
      )
    }
    
    
    entropy.balance.error = FALSE
    
    # perform entropy balancing
    # trim weights
    # then combine weights in to data frame
    # note weights and observations must be in the same order
    for(i in 1:arms){
      # do entropy balancing and get weights
      tryCatch({ 
        #code to try here:
        assign(paste("eb.out", i, sep=""),
               ebalance(Treatment = c(get(paste("txgroup", i, sep=""))$treatment, 
                                      rep(1, times = 10000)),
                        X = rbind((get(paste("txgroup", i, sep="")))[,1:4], 
                                  external_sample[,1:4]),
                        constraint.tolerance = 2
               )
        )}, error = function(e) {
          #error handling here
          entropy.balance.error <<- TRUE; 
          print(e)
      }
      )
    }
    print(entropy.balance.error)
    
    if(entropy.balance.error == FALSE){
      for(i in 1:arms){
        try(
          # trim weights
          assign(paste("eb.out", i, sep=""),
                 ebalance.trim(get(paste("eb.out", i, sep="")))
          )
        )
        # add weights to data
        assign(paste("txgroup", i, sep=""),
               (get(paste("txgroup", i, sep="")) %>% 
                  mutate(treatment = txNums[i],
                         weight = (get(paste("eb.out", i, sep="")))$w)
               )
        )
      }
    }
  
    
    if(entropy.balance.error == TRUE){
      base_weight = rep(1, NperInterim)
      for(i in 1:arms){
        assign(paste("txgroup", i, sep=""),
               (get(paste("txgroup", i, sep="")) %>% 
                  mutate(treatment = txNums[i],
                         weight = base_weight)
               )
        )
      }
    } 
    
    
    # create complete dataset with weights
    dataset = rbind(txgroup1, txgroup2)
    if(arms > 2){
      for(i in 3:arms){
        dataset = rbind(dataset, get(paste("txgroup", i, sep="")))
      }
    }
    
    
    ## 3.3) Fit Bayesian Model
    if(arms > 1){
    #select columns of treatment indicators for active treatments only
    x_model = as.numeric(dataset$treatment == txNums[1])
    for(i in 2:arms){
      x_model = cbind(x_model, as.numeric(dataset$treatment == txNums[i]))
    }
    
    #precision prior calc to reflect total number of participants in tx group so far
    print(thetapriors[, 2])
    prec_prior = 1/thetapriors[, 2]
    print(prec_prior)
    for(i in 1:length(thetapriors[, 2])){
      if(sampleSize[i] > NperInterim){
        prec_prior[i] = (sampleSize[i]/(2*NperInterim)) * prec_prior[i]
      }
    }
    print(prec_prior)
    
    # 3.3.1) write JAGS model for estimating posterior and comparing tx outcomes
    model_filename = paste("JAGSmodel_", trialName,   "_NperI", 
                           NperInterim, "_Arms", Arms, "_Thresh", supThresh, 
                           futThresh, ".txt", sep="")
    sink(model_filename)
    cat("
          model{
            for(i in 1:N){
              y[i] ~ dnorm(mu[i], prec * weight[i])
              mu[i] <- inprod(x[i,], theta[])
            }
            for(j in 1:arms){
              theta[j] ~ dnorm(mu_prior[j], prec_prior[j])
            }
            prec ~ dgamma(a, b)
          }
          ", fill = TRUE)
    sink()
    
    # 3.3.2) list of input data for model
    data1 = list(N = NperInterim*arms, 
                 x = x_model, 
                 y = dataset$outcome, 
                 weight = dataset$weight, 
                 arms = arms, 
                 mu_prior = thetapriors[, 1], 
                 prec_prior = prec_prior,
                 a = sigmaprior[1], 
                 b = sigmaprior[2]
    )
  
    # 3.3.3) initialize model hyper-parameters
    chains = 3
    inits1 = rep(list(list(theta = rep(0.5, arms), prec = 1)), chains)
    iter = 11000
    burnin = 1000
    thin = 2
    # parameters to track
    params = c("prec", "theta")
    
    
    # 3.3.4) fit JAGS model
    model = jags(data1, inits1, params, 
                 model_filename,
                 n.chains = chains, n.iter = iter, n.burnin = burnin, 
                 n.thin = thin, DIC = F)
    
    ## 3.4) Check Model Convergence
    modelConverge = checkModelConvergence(model, trialName, trialNum, k)
    
    # 3.4.1) If model did not converge:
    if(!modelConverge){
      print("first convergence check failed")
      
      # update model hyper-parameters
      iter = 32000
      burnin = 2000
      thin = 3
      
      # re-fit JAGS model
      model = jags(data1, inits1, params, 
                   model_filename,
                   n.chains = chains, n.iter = iter, n.burnin = burnin, 
                   n.thin = thin, DIC = F)
      
      # re-check model convergence
      modelConverge = checkModelConvergence(model, trialName, trialNum, k)
      
    }
    
    # 3.4.2) If model still did not converge:
    if(!modelConverge){
      # output autocorrelation and time series graphs to a pdf
      modelDiagnosticGraphs(model, trialName, trialNum, k)
    }
    
    ## 3.5) Analysis Summary
    
    # create analysis summary table (summary of mean outcomes)
    summary = trialSummaryNormal(model, txNames, arms, trialNum, txNums, k, 
                                 refFlag, theta, sampleSize, supThresh, 
                                 futThresh)
    summary[, 9] = as.numeric(summary[, 9]) + 10
    means = as.numeric(summary[, 6])
    vars = as.numeric(summary[, 7])^2
    superiority_probs = as.numeric(summary[, 8])
    cbind(summary, 
          `Entropy Error` = 
            as.numeric(rep(entropy.balance.error, arms))) -> summary
    
    # get summary of population variance estimate 
    varSummary = trialVarianceSummary(model)
    
    # add trial number and interim analysis number to summary
    varSummary = cbind(`Trial Number` = trialNum,
                       `Interim Analysis` = k,
                       varSummary)
    
    # 3.6) Save Summary Table to a CSV file
    write.table(summary, 
                file = paste("platformTrialSummary_", trialName, 
                             "_NperI", NperInterim, "_Arms", Arms, "_Thresh", 
                             supThresh, futThresh, 
                             ".csv", sep = ""), 
                row.names = FALSE, col.names = (k == 1 && trialNum == 1), na = "", 
                sep = ",", append = (k > 1 | trialNum > 1))
    
    # save variance summary
    write.table(varSummary, 
                file = paste("platformTrialVaraince_", trialName, 
                             "_NperI", NperInterim, "_Arms", Arms, "_Thresh", 
                             supThresh, futThresh,
                             ".csv", sep = ""), 
                row.names = FALSE, col.names = (k == 1 && trialNum == 1), na = "", 
                sep = ",", append = (k > 1 | trialNum > 1))
    
    # save distribution of outcomes
    tosave = as.matrix(t(
      c(trialNum, k, length(dataset$outcome), mean(dataset$outcome), 
        quantile(dataset$outcome, probs = c(0, 0.25, 0.5, 0.75, 1)))))
    
    colnames(tosave) = c("Trial.Number", "Analysis.Number", "N", "Mean.Outcome",
                         "Min", "25%", "Median", "75%", "Max")
    
    write.table(tosave, 
                file = paste("outcomes_", trialName, 
                             "_NperI", NperInterim, "_Arms", Arms, "_Thresh", 
                             supThresh, futThresh,
                             ".csv", sep = ""), 
                row.names = FALSE, col.names = (k == 1 && trialNum == 1), na = "",
                sep = ",", append = (k > 1 | trialNum > 1))
      
    ## 3.7) Update Priors
    # update priors for next iteration based on current data
    #   - update estimated mean and variance for theta priors (normal distribution)
    #   - update shape parameters a, b for sigma prior (gamma distribution)
    #
    # mean     
    thetapriors[, 1] = means
    # variance
    thetapriors[, 2] = vars
    
    # solve for gamma shape parameters, a and b
    # mean = a / b
    # variance = a / b ^2
    # >> a = mean^2/variance
    # >> b = mean/variance
    sigmamean = as.numeric(varSummary[1, 5])
    sigmavar = as.numeric(varSummary[1, 6])
    # a
    sigmaprior[1] = sigmamean^2/sigmavar
    # b
    sigmaprior[2] = sigmamean/sigmavar
    
    
    ## 3.8) Drop Treatment(s) for Futility
    #   Check futility criteria, if met, remove this treatment arm
    #   1) remove from theta vector
    #   2) remove from thetapriors vector
    #   3) remove from txNames vector and txNums vector
    #   4) update number of treatment arms
    #   5) update sample size tracker
    futile = which(as.numeric(superiority_probs) < futThresh)
    if(length(futile) > 0){
      theta = theta[-futile]
      thetapriors = thetapriors[-futile, ]
      txNames = txNames[-futile]
      txNums = txNums[-futile]
      means = means[-futile]
      superiority_probs = superiority_probs[-futile]
      arms = arms - length(futile)
      sampleSize = sampleSize[-futile]
    }
    
    ## 3.9) Superiority Switch
    # check superiority criteria
    #   - if met, move this treatment to position 1 so it will become the 
    #     new reference group
    #   - if more than one tx meet criteria, one with highest estimated
    #     success probability is used as new reference
    #   1) change order of theta vector
    #   2) change order of thetapriors vector
    #   3) change order of txNames vector and TxNums vector 
    #       (add and remove '(Ref)' label for names)
    #   4) change order of sample size tracker
    superior = which(as.numeric(superiority_probs) > supThresh)
    if(length(superior) > 1){
      superior = superior[which.max(means[superior])]
    }
    if(length(superior) == 1){
      theta = c(theta[superior], theta[-superior])
      thetapriors = rbind(thetapriors[superior, ], thetapriors[-superior, ])
      txNames = c(paste(txNames[superior], "(Ref)", sep = " ")
                  , txNames[-superior])
      txNames[2] = substr(txNames[2], 1, str_length(txNames[2])-6)
      txNums = c(txNums[superior], txNums[-superior])
      sampleSize = c(sampleSize[superior], sampleSize[-superior])
    }
    
    
    ## 3.10) Add in New Treatment(s)
    #   - add more if any treatments have been dropped for futility
    #   1) add to theta vector
    #   2) add to thetapriors vector
    #   3) add to txNames vector and txNums vector
    #   4) update number of treatment arms
    #   5) update number of adaptations tested
    if(length(futile) > 0){
      #for each adaptation dropped for futility
      for(i in 1:length(futile)){
        # if have not reached max number of adaptations
        if (numAdapts + 1 <= length(txEffects)){
          # add new adaptation with next value from possible theta parameters
          theta = c(theta, txEffects[numAdapts + 1])
          thetapriors = rbind(thetapriors, c(0, 10))  
          numAdapts = numAdapts + 1
          arms = arms + 1
          txNames = c(txNames, paste("Adaptation", numAdapts, sep = " "))
          txNums = c(txNums, numAdapts)
          sampleSize = c(sampleSize, 0)
        }
      }
    }
    
    refFlag = c(1, rep(0, arms-1))
    }
      
    if(arms == 1){
      
      #save summary table
      summary = as.matrix(t(c(trialNum, k, txNums, 1, txNames, 
                  round(sum(dataset$outcome*dataset$weight)/sum(dataset$weight), 3),
                  NA, NA, theta, NA, sampleSize, NA, NA, 
                  as.numeric(entropy.balance.error))))
      write.table(summary, 
                  file = paste("platformTrialSummary_", trialName, 
                               "_NperI", NperInterim, "_Arms", Arms, "_Thresh", 
                               supThresh, futThresh, 
                               ".csv", sep = ""), 
                  row.names = FALSE, col.names = FALSE, na = "", 
                  sep = ",", append = TRUE)
      
      # save distribution of outcomes
      tosave = as.matrix(t(c(trialNum, k, length(dataset$outcome), 
                            mean(dataset$outcome), 
                            quantile(dataset$outcome, 
                                     probs = c(0, 0.25, 0.5, 0.75, 1)))))
      
      write.table(tosave, 
                  file = paste("outcomes_", trialName, 
                               "_NperI", NperInterim, "_Arms", Arms, "_Thresh", 
                               supThresh, futThresh,
                               ".csv", sep = ""), 
                  row.names = FALSE, col.names = FALSE, na = "",
                  sep = ",", append = TRUE)
    }
    
    ## 3.11) Update Counters
    
    # once last adaptation is added, start incrementing counter
    # want to allow up to 3 interim analyses after last adaptation is added
    if(numAdapts == length(txEffects)){
      maxIteration = maxIteration + 1
    }
    
    # update total sample size
    totalSample = totalSample + (arms * NperInterim)
    
    # update counter for number of interim analyses
    k = k + 1
  }
  
  ## 4) Record Final Sample Size
  
  # since total sample size is updated at the end of the loop in anticipation 
  # of the next iteration, subtract last update in order to get true ending
  # sample size
  totalSample = totalSample - (arms * NperInterim)
  totalSample = cbind(trial = trialNum, sampleSize = totalSample)
  
  # record final sample size, save to csv file
  write.table(totalSample, sep = ",", append = (trialNum > 1), 
              row.names = FALSE, col.names = (trialNum == 1), 
              file = paste("totalSampleSizes_", trialName,   
                           "_NperI", NperInterim, "_Arms", Arms, "_Thresh", 
                           supThresh, futThresh,
                           ".csv", sep = ""))
  
  # IF randomly generated order of txEffects, record the order
  if(str_detect(trialName, "Random")){
    write.table(t(txEffects), sep = ",", append = (trialNum > 1), 
                row.names = FALSE, col.names = FALSE, 
                file = paste("txEffectOrder_", trialName,   
                             "_NperI", NperInterim, "_Arms", Arms, "_Thresh", 
                             supThresh, futThresh,
                             ".csv", sep = ""))
  }
}

################################################################################

