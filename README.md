# HRMC
High-rank matrix completion for gene prioritisation

Matlab scripts to run the algorithm.

* HRMC.m contains the implementation of our algorithm.
* example_script.m contains an example script of how to run the complete model.

The following data in matlab file-format for using this code can be found at http://www.paccanarolab.org//hrmc-gene/

``` Matlab
% DGAM2017: gene-disease association 2017.
% DSIM2017: Caniza et. al. semantic similarity 2017.
% PPI: HPRD protein interaction network.
% genes_names: entrez IDs.
% disease_MIM: IDs.
```
The usage of the HRMC model is very simple, explained step-by-step below:

## Step 1: Set algorithm parameters
``` Matlab
gamma = 10^4;
variance = 0.01;
tolX = 1e-2;
maxiter = 100;
```

## Step 2: Column model HRMC-c

It uses the graph regularization on the semantic similarities between diseases.

``` Matlab
C = HRMC(DGAM2017, DSIM2017,...
            0.5,...
            1, 0.5,...
           gamma, variance,...
           tolX, maxiter);
       
HRMCc = DGAM2017 * C;
```
## Step 3: Row model HRMC-r

It uses the graph regularization on the human PPI network.

``` Matlab
R = HRMC(DGAM2017', PPI,...
            0.5,...
           0.5, 0.5,...
            gamma, variance,...
           tolX, maxiter);

HRMCr = (DGAM2017' * R)';
```
## Step 4: Final linear combination
``` Matlab
p = 0.7;
Xhat = p*HRMCc + (1-p)*HRMCr; 
``` Matlab
