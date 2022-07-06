OFBM tools
===

## Description
OFBM tools is a multivariate sel-similarity analysis matlab package. It permits to estimate and group the values of scaling exponents across multivariate time series modeled by an operator fractional Brownian motion (ofBm). The estimation is based on linear regressions of wavelet eigenvalues across scales. The clustering strategy relies on successive pairwise tests between ordered estimates. Two methods are proposed to estimate the paramters of the pairwise tests.

## Requirements
The following matlab package is required: [WLBMF](https://www.irit.fr/~Herwig.Wendt/software.html).

## References
  - [Lucas et al., 2021](https://eurasip.org/Proceedings/Eusipco/Eusipco2021/pdfs/0001960.pdf)
  - [Lucas et al., 2022](https://ieeexplore.ieee.org/document/9747448)
  - 
  
## Quick start
The basic syntax to run OFBM Tools is as follows:

```
% parameters of the estimation
paramsEst.Nwt = 2 ; paramsEst.FigNum = 10 ; paramsEst.wtype = 1 ;
paramsEst.j1 = 8; paramsEst.j2 = 11; paramsEst.Jref = paramsEst.j2 ; 
paramsEst.NB = 0; paramsEst.LB = 0;
% return self-similarity exponent estimation
[est,estbc] = OFBM_estimBC_BS(data,paramsEst) ;
```
The parameters to take into account in the input structure `paramsEst` of OFBM_estimBC_BS are:
  - `j1`, first scale for analysis;
  - `j2`, last scale for analysis;
  - `Jref`, reference scale under which several wavelet spectra are computed with the same number wavelet coeficients;
  - `wtype`, kind of weighting used in the linear regressions:
    - 0 for no weighting  (i.e., uniform weights);
    - 1 for 1/nj weights  (suitable for fully Gaussian data);
    - 2 to use variance of the estimates;
  - `NB`, number of bootstrap resampling;
  - `LB`, number of blocks for the bootstrap resampling;
  - `Nwt`, number of vanishing moments of the wavelet;
  - `FBM`, type of the process:
    - 1 for operator fractional Brownian motion (ofBm);
    - 0 for operator fractional Gaussian noise (ofGn).

The main parameters contained in the structures `est` and `estbc` returned by OFBM_estimBC_BS are:
  - `est.hU`, matrix of univariate-like estimates of the self similarity exponents;
  - `est.h`, classical multivariate estimates of the self similarity exponents;
  - `estbc.h`, bias corrected multivariate estimates of the self similarity exponents.
    
The clustering of the scaling exponents need to run the estimation with adapted parameters:
```
% parameters of the estimation for clustering 
% bootstrap estimates are needed for the pairwise tests
paramsEst.NB = 500; paramsEst.LB = 2*params.Nwt; 
% return self-similarity exponent estimation
[est,estbc] = OFBM_estimBC_BS(data,paramsEst) ;

% testing procedure
alpha = 0.05; % significance level
paramsTest.P = size(data,1); paramsTest.NB = paramsEst.NB;
estT = OFBM_estimBC_BS_test(estbc,alpha,paramsTest);
% cluster the self-similarity exponents
[nbcluster,cluster] = successiveTestClustering(estT.decsortHocpw);
```
The parameters to take into account in the input structure `paramsTest` of OFBM_estimBC_BS_test are:
  - `P`, the number of components;
  - `NB`, number of bootstrap resampling.


Another clustering method is available:
```
% cluster the self-similarity exponents
[nbcluster,cluster] = successiveTestClustering(estT.decsortHocpw_v2);
```

Examples with simulated ofBm can also be found in the `examples` folder.
