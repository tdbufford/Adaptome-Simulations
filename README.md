# Adaptome-Simulations

The code here contains simulations to demonstrate feasibility and performance of a statistical methodology, called the "Adaptome." This novel statisical method combines existing methods used in platform clinical trials for comparison of multiple treatments via interim analyses with existing methods used in analyzing obervational data to balance covariate distributions among groups. The chosen covariate balaincing method is called "entropy balancing."

The purpose of this method is to use data collected in real-world settings to improve health outcomes by identifying beneficial adaptations to health interventions that arise during intervention implementation.

The files with names that begin "PlatformTrial_" define functions that run simulated versions of a platform trial. Data are generated for mulitple intervention groups, intervention effectiveness is assessed, inferior intervention adaptations are dropped and replaced with new adaptations, and then the process repeats. The trial ends after all possible intervention adaptations have been assessed. Simulation results are output to .csv files.

The suffix "Entropy" in the file name indicates that the method used to assess intervention effectiveness includes entropy balancing.
The suffix "EntropyCovars" in the file name indicates that the method used to assess intervention effectiveness includes entropy balancing and also adjusts for covariates in the Bayesian model.
the suffix "noEntropy" in the file name indicates that the method used to assess intervention effectiveness does not include entropy balancig, or any adjustment for covariates.

The files  with names that begin "helper_" contain helper functions called by the platform trial functions.

The files with names that begin "RUN" contain code that runs the simulated platform trial by calling the functions defined in the other files. Various versions of the simulation are run by varying the inputs. 

The files with names that begin "SimResults" contain code to read the .csv files output from the platform trial functions and compute performance metrics.

The files with the suffix "SchoolData" contain versions of the simulation where the generated data structure is designed to mimick data collected from a health intervention implemented in a large school district. Entropy balancing is used in these simulated platform trials. This version of the simulation demonstrates possible improvements to student outcomes when data are analyzed intermittently using the Adaptome statistical method. 
