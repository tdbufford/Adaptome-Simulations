##### checkModelConvergence function #############################################
# check model convergence using:
#     1) Gelman-Ruben estimated potential scale reduction factor (shrink factor)
#     2) Autocorrelation assessment 
#     3) Cross-correlation assessment (between parameters)
#     diagnostics are calculated for first MCMC chain
# @param model - output from running the JAGS model
checkModelConvergence = function(model, trialName, trialNum, k){
  #get sims array from model - this is samples from the posterior for each param
  #the dimensions are: (iterations, chains, parameters)
  array = model$BUGSoutput$sims.array
  
  #get numbermodel parameters to track (number treatment arms + number covariates)
  nparams = dim(array)[3]
  
  #set convergence check flag - if any test does not pass, flag set to FALSE
  converge = TRUE
  
  #set warning message
  warning = paste(trialName, " trial ", trialNum, ", Interim analysis ", k, 
                  ". Model convergence check failed", sep = "")
  
  #check autocorrelation for chain 1 for each parameter
  #set convergence flag to FALSE if autocorrelation exceeds threshold for any parameter
  auto = apply(array[, 1, ], 2, function(x) acf(x, plot = FALSE)$acf)
  if(sum(auto[2,] > 0.4) > 0){converge = FALSE; 
  warning = paste(warning, ", model has high autocorrelation",
                  sep = "")}
  if(sum(auto[5,] > 0.1) > 0){converge = FALSE; 
  warning = paste(warning, ", model has high autocorrelation",
                  sep = "")}
  # if(sum(auto[10,] > 0.05) >0){converge = FALSE; 
  #                 warning = paste(warning, ", model has high autocorrelation",
  #                 sep = "")}
  
  #check cross-correlation (correlation between parameters)
  #set convergence flag to FALSE if any two parameters have correlation greater than 0.5
  corr = cor(array[, 1, ], method = "pearson")
  if(sum(abs(corr) > 0.5) > nparams){converge = FALSE; 
  warning = paste(warning, ", model has cross correlation",
                  sep = "")}
  
  #calculate the shrink factor
  #As the simulation converges, the shrink factor declines to 1
  shrink = gelman.diag(model$BUGSoutput, confidence = 0.95, transform = FALSE,
                       autoburnin = FALSE, multivariate = TRUE)
  if(sum(shrink$psrf > 1.1) > 0){converge = FALSE; 
  warning = paste(warning, ", model has high shrinkage factor", 
                  sep = "")}
  if(shrink$mpsrf > 1.1){converge = FALSE; 
  warning = paste(warning, ", model has high shrinkage factor",
                  sep = "")}
  
  
  if(!converge){
    filename = paste("ModelConvergence_", trialName, ".txt", sep = "")
    sink(file = filename, append = TRUE)
    cat(paste(warning, ". ", Sys.time(), sep = ""), append = TRUE, sep = "\n")
    cat("\n", append = TRUE, sep = "\n")
    sink()
  }
  
  return(converge)
}