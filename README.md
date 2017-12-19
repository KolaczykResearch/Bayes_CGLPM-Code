% This contains a data and functions folder with scripts to implement the 
% experiment in Ahelegbey D., Carvalho L., Kolaczyk E. (2017). 
% "A Bayesian Covariance Graphical And Latent Position Model For 
% Multivariate Financial Time Series"
% =======================================================================
%
% Demo.m - a scrip to implement the simulation experiment in in Section 3.
% 
% Generate_Data.m - script used to generate the data for the simulation experiment in Section 3. 
% 
% All data files and script file for data generation are available in the data folder.
% 
% BCGLPM_MCMC.m - scrip for MCMC sampling of the Bayesian covariance graphical structure (G) and the 
% latent positions (U) of the nodes.
% 
% COMPUTE_SS_MATRIX.m - scrip for computing the unnormalized covariance matrix of Y conditional on X.
% 
% GRID_SEARCH.m - scrip for the grid search of the ridge parameter
% 
% Matrix_Bingham_VMF_Gibbs.m - scrip for Gibbs sampling of the latent positions (U) from a matrix 
% Bingham-von Mises-Fisher distribution.
% 
% Sample_Z_FC.m - scrip for sampling the latent Z following Hoff's rstiefel R package
% 
% PSRF_CONVERGENCE.m - scrip for computing the Gelman and Rubin (1992) potential scale reduction factor (PSRF) 
% convergence diagnostics of the MCMC.
% 
% Procrustean_Analysis.m - scrip for Procrustes transformation of U_1 to U_0 (with U_0 as the target)
% 
% gigrnd.m - scrip for implementing the Devroye (2014) algorithm for sampling from the generalized 
% inverse Gaussian (GIG) distribution
% 
% All these files are contained in the function folder.
