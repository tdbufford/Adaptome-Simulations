####### trialVarianceSummary function #####################################################
# Get posterior estimate for data variance
# and create summary table for interim analysis 
# this summary is only valid for data with continuous outcomes
# (otherwise variance is not estimated)
# 
# @param model - output from running the JAGS model
# @return summary - table of values:
#                   posterior mean, sd, and superiority prob 
#                   one row for each parameter
trialVarianceSummary = function(model){
  #get sims array from model - this is samples from the posterior for each param
  #the dimensions are: (iterations, chains, parameters)
  array = model$BUGSoutput$sims.array
  
  #get number of chains
  nchains = dim(array)[2]
  
  # #get number of parameters (last one is the precision) 
  # nparams = dim(array)[3]
  
  #select only posterior samples for precision parameter
  #array = array[, , nparams]
  array = array[, , 1]
  
  #reshape model output
  #output: rows = iterations (chains combined)
  output = as.matrix(array[, 1])
  for (i in 2:nchains){
    output = rbind(output, as.matrix(array[, i]))
  }
  
  #get total number of samples from posterior
  nsamples = nrow(output)
  
  #calc posterior mean
  yvar = 1/output
  posterior_mean = round(mean(yvar), 3)

  #calc posterior sd
  posterior_sd = round(sd(yvar), 3)
  
  prec_mean = mean(output)
  prec_var = var(output)
  
  summary = cbind(posterior_mean, posterior_sd, prec_mean, prec_var)
  
  #add column names
  colnames(summary) = c("Population Variance Estimate (Posterior Mean)", 
                        "SD of Estimate", "Precision_mean", "Precision_var")
  
  return(summary)
}