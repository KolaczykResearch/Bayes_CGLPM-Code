% This code implements the Simulation Experiment in 
% Ahelegbey D., Carvalho L., Kolaczyk E. (2017). 
% "A Bayesian Covariance Graphical And Latent Position Model
% For Multivariate Financial Time Series"
%========================================================================
close all; clc; clear;

% Add path of functions, data and figures
addpath('functions');    addpath('data');    

%----- GENERATE SIMULATED DATA
n = 50;             % number of variables (multiple of 50 e.g 50, 100, 150)
T = 2*n;            % number of observations (multiple of n e.g 2n, 10n)
lag = 1;            % lag order (0, 1)
[Endo, Exo, SigTrue]  = Generate_Data(n, T, lag);
AdjTrue = double(abs(SigTrue)>1e-4) - eye(n);

%------ PRELIMINARIES
nsimu  = 1e4;       % number of Simulation
v_0 = 0.02;         % Hyper-parameter v_0
h  = 1/v_0;         % Hyper-parameter h

%----- SAMPLING NETWORK AND LATENT POSITIONS 
lag0 = 0;
lag1 = 1;
M_0 = BCGLPM_MCMC(Endo, Exo, lag0, nsimu, h, v_0);
M_1 = BCGLPM_MCMC(Endo, Exo, lag1, nsimu, h, v_0);

%----- EVALUATE PROCRUSTEAN DISTANCE
Proc = Procrustean_Analysis(M_0.U, M_1.U);