
Bayesian Covariance Graphical And Latent Position Model
===============================

Matlab scripts to implement the experiment in the paper "A Bayesian Covariance Graphical And Latent Position Model For Multivariate Financial Time Series" by Ahelegbey D., Carvalho L., Kolaczyk E.

  Demo.m 
  - scrip to implement the simulation experiment in Section 3.

data folder
--
  Generate_Data.m 
  - script used to generate the data for the simulation experiment in Section 3. 


functions folder
--

  BCGLPM_MCMC.m
  - scrip to implement MCMC sampling of the covariance graphical structure (G) and the latent positions (U) of the nodes.

  COMPUTE_SS_MATRIX.m 
  - scrip for computing the unnormalized covariance matrix of Y conditional on X.

  GRID_SEARCH.m 
  - scrip for the grid search of the ridge parameter

  Matrix_Bingham_VMF_Gibbs.m 
  - scrip to implement Gibbs sampling of the latent positions (U) from a matrix Bingham-von Mises-Fisher distribution.

  Sample_Z_FC.m 
  - scrip for sampling the latent Z following Hoff's rstiefel R package

  PSRF_CONVERGENCE.m 
  - scrip for computing the Gelman and Rubin (1992) potential scale reduction factor (PSRF) convergence diagnostics of the MCMC.

  gigrnd.m 
  - scrip to implement the Devroye (2014) algorithm for sampling from the generalized inverse Gaussian (GIG) distribution
