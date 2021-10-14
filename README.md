OFBM tools
===

## Description
OFBM tools is a multivariate sel-similarity analysis matlab package. It permits to estimate and group the values of scaling exponents across multivariate time series modeled by an operator fractional Brownian motion (ofBm). The estimation is based on linear regressions of wavelet eigenvalues across scales. The clustering strategy relies on successive pairwise tests between close estimates.

## Requirements
The following matlab package is required: [WLBMF](https://www.irit.fr/~Herwig.Wendt/software.html).

## References
  - [Lucas et al., 2021](https://www.irit.fr/~Herwig.Wendt/data/LucasEUSIPCO2021.pdf)
  
## Quick start
The basic syntax to run OFBM Tools is as follows:

```
% parameters of the estimation
paramsEst = params; % params contains the parameters of the ofBm generated in the `data` folder
paramsEst.Nwt = 2 ; paramsEst.FigNum = 10 ; paramsEst.wtype = 1 ;
paramsEst.j1 = 8; paramsEst.j2 = 11; paramsEst.Jref = paramsEst.j2 ; 
paramsEst.NB = 0; paramsEst.LB = 0;
% return self-similarity exponent estimation
[est,estbc] = OFBM_estimBC_BS(data,paramsEst) ;
```

The clustering of the scaling exponents need to run the estimation with adapted parameters:
```
% parameters of the estimation for clustering (bootstrap estimates are needed for the pairwise tests)
paramsEst.NB = 500; paramsEst.LB = 2*params.Nwt; 
% return self-similarity exponent estimation
[est,estbc] = OFBM_estimBC_BS(data,paramsEst) ;

% testing procedure
alpha = 0.05;
estT = OFBM_estimBC_BS_test(estbc,alpha,params);
% clustering of the self-similarity exponents
[nbcluster, cluster] = successiveTestClustering(estT.decsortHocpw);
```

The main parameters to take into account in the structure `params` are:

  - `R`, the number of realizations of the Monte Carlo vector;
  - `Jref`: reference scale under which several wavelet spectra are computed with the same number wavelet coeficients
  - `j2`: last scale for analysis
  - `j1`: first scale for analysis
  - `wtype`, kind of weighting used in the linear regressions:
    - no weigthing  (ie uniform weights)
    - 1/nj weights  (suitable for fully Gaussian data)
    - use variance of the estimates
  - `NB`: number of bootstrap resampling
  - `LB`: number of blocks for the bootstrap resampling
  
Examples with simulated ofBm can also be found in the `examples` folder.
