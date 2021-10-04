OFBM tools
===

## Description
OFBM tools is a multivairate sel-similarity analysis matlab package. It permits to estimate and group the values of scaling exponents across multivariate time series modeled by an operator fractional Brownian motion (ofBm).

## Requirements
The following matlab package is required: [WLBMF](https://www.irit.fr/~Herwig.Wendt/software.html).

## References
  - [Lucas et al., 2021](https://www.irit.fr/~Herwig.Wendt/data/LucasEUSIPCO2021.pdf)
  
## Quick start
The basic syntax to estimate self-similarity exponents is as follows:

```
% parameters of the estimation
paramsEst = params; % params contains the parameters of the ofBm generated in the `data` folder
paramsEst.Nwt = 2 ; paramsEst.j1 = 5; paramsEst.j2 = 10; paramsEst.Jref = paramsEst.j2 ;
paramsEst.FigNum = 10 ; paramsEst.wtype = 1 ; paramsEst.NB = 0; paramsEst.LB = 0;
% return self-similarity exponent estimation
[est,estbc] = OFBM_estimBC_BS(data,paramsEst) ;
```

The basic syntax to cluster self-similarity exponents is as follows:
```
% parameters of the estimation
paramsEst = params; 
paramsEst.Nwt = 2 ; paramsEst.j1 = 8; paramsEst.j2 = 11; paramsEst.Jref = paramsEst.j2 ;
paramsEst.FigNum = 10 ; paramsEst.wtype = 1 ; 
paramsEst.NB = 500; paramsEst.LB = params.Nwt; % bootstrap estimates are needed for the pairwise tests
% return self-similarity exponent estimation
[est,estbc] = OFBM_estimBC_BS(data,paramsEst) ;

% Testing procedure
alpha = 0.05;
estT = OFBM_estimBC_BS_test(estbc,alpha,params);
% Clustering of the self-similarity exponents
[nbcluster, cluster] = successiveTestClustering(estT.decsortHocpw);
```

The main parameters to take into account in the structure `params` are:

  - `R`, the number of realizations of the Monte Carlo vector;
  - `Jref`: reference scale under which several wavelet spectra are computed with the same number wavelet coeficients
  - `j2`: last scale for analysis
  - `j1`: first scale for analysis
  - `wtype`, the kind of weighting used in the linear regressions:
    - no weigthing  (ie uniform weights)
    - 1/nj weights  (suitable for fully Gaussian data)
    - use variance of the estimates
  - `NB`: number of bootstrap resampling
  - `LB`: number of blocks for the bootstrap resampling
  
Examples with simulated ofBm can also be found in the `examples` folder.
