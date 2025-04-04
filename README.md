#################################################################################

Source code and data for the manuscript 
"Robust Bayesian Inference for Accelerated Failure Time Models with Skewed and Heavy-Tailed Survival Data"
by Narjes Gilani, Mehrdad Naderi, and Reza Pourmousa
#################################################################################

# Author responsible for the code #
For questions, comments, or remarks about the code please contact the responsible author, Narjes Gilani (Narges_Gilani@math.uk.ac.ir).

# Configurations #
Please copy the files to the "current working directory" of the R.

# Descriptions of the codes #
 1. MH_nu.r: Sampling from the posterior distribution of nu under different prior specifications.
 2. NMV_CLR_Bayes.r:  A parametric fit for censored linear regression models based on NMVM distributions, from a Bayesian perspective
 3. Gibbs-sampling.r: Sampling from the posterior distribution of various parameters.
 4. NMVMdist: PDF and CDF of the NMVM family of distributions
 5. Example code to use the program.
 
# Notes #
  1. The "ggplot2","gridExtra","truncnorm","ghyp","cli" and "tcltk" R packages are required;
  2. The "Breast cancer" data set is available in the R package "Survival";
  3. The "Ambulatory expenditures" data set is available in the R package "ssmrob";
# parameters

#'@param{
#' y~ denotes the response vector, which may be right-censored or left-censored

#' x~ A matrix or vector of explanatory (predictor) variables.

#' Family~  Distribution to be used: "Normal", "VG","GHST", "NMVBS", "H", "NMVL", "AL".

#' cens~ "Left" for left-censoring, "Right" for right-censoring.

#' cc~ Censoring indicator vector: 1= observed, 0= censored.

#' influence~ "TRUE" or "FALSE". Indicates if the divergence measures (KL divergence, J, L and Chi Distance) should be computed.

#' spacing~ The 'spacing' parameter should only be set when either 'influence' or 'criteria' is TRUE, defining the lag between retained post-burn-in/thinning samples used for computation (spacing=1 uses all samples)."

#' criteria~ "TRUE" or "FALSE". Indicates if model selection criteria (LPML, DIC, EAIC, EBIC and WAIC) should be computed.

#' influence~ "TRUE" or "FALSE". Indicates if the divergence measures.

#' (KL divergence, J, L and Chi Distance) should be computed.

#' hyper_set~ enables assigning new values to the model's hyperparameters.

#' prior~ Possible prior distributions that may be considered for parameter nu include:  
"Exp" (exponential distribution), "Unif"(Uniforme distribution), and "Hierar_1" (Hierarchical prior distribution).
#' hyper~ Value of hyperparameter for the exponential prior.

#' n.thin~ Lag for posterior sample thta reduces autocorrelation between MCMC samples by retaining every nth iteration.  

#' burnin~ Burn-in phase in Bayesian MCMC discards initial samples to ensure convergence before analysis.

#' n.iter~ The number of iterations to be considered (before burnin and thinning).

#' n.chains~ The number of chains to be considered.

#' sample.out~ If "TRUE", all the posterior chains are stored for posterior analysis.

#' level~  indicates the percentage of censorship.
