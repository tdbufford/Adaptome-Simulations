##### modelDiagnosticGraphs function #############################################
# creates plots to check model convergence
#     includes autocorrelation plots and time series plots 
#     for first MCMC chain
# plots are output to a pdf file defined in simulation set up
#     "platformTrialModelConvergence.pdf"
# @param model - output from running the JAGS model
modelDiagnosticGraphs = function(model, trialName, trialNum, k){
  # set file for saving plots 
  filename = paste("ModelConvergence_", trialName, trialNum, "_", k, ".pdf", sep = "" )
  pdf(file = filename)
  
  #get sims array from model - this is samples from the posterior for each param
  #the dimensions are: (iterations, chains, parameters)
  array = model$BUGSoutput$sims.array
  
  #get number of treatment arms (number of model parameters to track)
  narms = dim(array)[3]
  
  #Plot autocorrelation for chain 1 for each parameter
  par(mfrow = c(ceiling(narms/2), 2))
  for(i in 1:narms){
    acf(array[, 1, i], main = paste("Analysis", k, "- Param", i, sep = " ")) 
  }
  #Plot time series for chain 1 for each parameter
  par(mfrow = c(ceiling(narms/2), 2))
  for(i in 1:narms){
    plot(1:1000, array[1:1000, 1, i], type="l", xlab="iteration", ylab="",
         main = paste("Analysis", k, "- Param", i, sep = " "))
    # lines(4000:5000, array[4000:5000, 2, 4], type="l", xlab="iteration", ylab="", col="blue")
    # lines(4000:5000, array[4000:5000, 3, 4], type="l", xlab="iteration", ylab="", col="magenta")
    # lines(4000:5000, array[4000:5000, 4, 4], type="l", xlab="iteration", ylab="", col="purple")
  }
  chart.Correlation(array[, 1, ], histogram=TRUE, pch=12)
  #close file
  dev.off()
}