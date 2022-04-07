####### trialSummaryNormal function #####################################################
# calculate prob for superiority, P(theta[i]>theta[1] | data)
# and create summary table of interim analysis 
# including posterior means and sds of estimated thetas
# 
# @param model - output from running the JAGS model
# @param names - vector of treatment names
# @param arms - number of active treatment arms in trial
# @return summary - table of values:
#                   posterior mean, sd, and superiority prob 
#                   one row for each parameter
trialSummaryNormal = function(model, txNames, arms, trialNum, txNums, k, refFlag,
                              theta, sampleSize, supThresh, futThresh){
  #get sims array from model - this is samples from the posterior for each param
  #the dimensions are: (iterations, chains, parameters)
  array = model$BUGSoutput$sims.array
  
  #get number of chains
  nchains = dim(array)[2]
  
  #get number of treatment arms 
  #arms = dim(array)[3]
  
  array = array[, , 2:(arms+1)]
  
  if(length(txNames) != arms){print("List of names does not match number of 
                                   treatment arms")}
  
  #reshape model output
  #output: rows = iterations (chains combined), cols = theta
  output = array[, 1, ]
  for (i in 2:nchains){
    output = rbind(output, array[, i, ])
  }
  
  #get total number of samples from posterior
  nsamples = nrow(output)
  
  #calc posterior mean
  posterior_means = apply(output, 2, mean)
  posterior_means = round(posterior_means, 3)
  
  #calc posterior sd
  posterior_sd = apply(output, 2, sd)
  posterior_sd = round(posterior_sd, 3)
  
  #calc probability of superiority compared to reference treatment
  #first treatment is assumed to be reference group
  prob_superior = rep(NA, arms)
  for (i in 2:arms){
    prob_superior[i] = round(length(which(output[,i]>output[,1]))/nsamples, 3)
  }
  
  summary = cbind(`Trial Number` = rep(trialNum, arms),
                  `Interim Analysis` = rep(k, arms), 
                  `Intervention Number` = txNums, 
                  `Reference Group Flag` = refFlag, 
                  `Intervention Name` = txNames, 
                  `Estimate` = posterior_means, 
                  `SD of Estimate` = posterior_sd, 
                  `Superiority Prob` = prob_superior,
                  `True Value` = theta, 
                  `Effect Size` = theta - theta[1], 
                  `Sample Size` = sampleSize,
                  `Superiority Criteria Met` = ifelse(prob_superior > supThresh, 1, 0),
                  `Inferiority Criteria Met` = ifelse(prob_superior < futThresh, 1, 0))
  
  # effect size should not be calculated for reference group
  summary[1, 10] = NA
  
  return(summary)
}