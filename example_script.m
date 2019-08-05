%% HRMC - High-rank matrix completion for gene prioritization
% Authors: Cheng Ye, Diego Galeano, Alberto Paccanaro.
% Copyright @ 2019. Code: Diego Galeano.

%% Load the data 
% Data is here: http://www.paccanarolab.org//hrmc-gene/
% DGAM2017: gene-disease association 2017.
% DSIM2017: Caniza et. al. semantic similarity 2017.
% PPI: HPRD protein interaction network.
% genes_names: entrez IDs.
% disease_MIM: IDs.
clear all; clc;
load('data');

%% Set model parameters
gamma = 10^4;
variance = 0.01;
tolX = 1e-2;
maxiter = 100;

%% Column model
tic;
C = HRMC(DGAM2017, DSIM2017,...
            0.5,...
            1, 0.5,...
           gamma, variance,...
           tolX, maxiter);
       
HRMCc = DGAM2017 * C;
toc;
%% Row model
tic;
R = HRMC(DGAM2017', PPI,...
            0.5,...
           0.5, 0.5,...
            gamma, variance,...
           tolX, maxiter);

HRMCr = (DGAM2017' * R)';
toc;

%% Linear combination 
p = 0.7;
Xhat = p*HRMCc + (1-p)*HRMCr; 


%% --------------   


