OFBM tools
===

## Description
OFBM tools is a multivariate self-similarity analysis matlab package. It permits to estimate and count the values of scaling exponents across multivariate time series modeled by an operator fractional Brownian motion (ofBm). The estimation is based on linear regressions of wavelet eigenvalues across scales. The clustering strategy relies on successive parametric pairwise hypothesis tests between ordered estimates: two bootstrap-based estimation methods are proposed for the test parameters.

## Requirements
The following matlab package is required: [WLBMF](https://www.irit.fr/~Herwig.Wendt/software.html).

## References
  - [Lucas et al., 2022](https://ieeexplore.ieee.org/document/9747448)
  - [Lucas et al., 2021](https://eurasip.org/Proceedings/Eusipco/Eusipco2021/pdfs/0001960.pdf)
  
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
  - `FBM`, type of the process:
    - 1 for operator fractional Brownian motion (ofBm);
    - 0 for operator fractional Gaussian noise (ofGn);
  - `Nwt`, number of vanishing moments of the wavelet;
  - `j1`, first scale for analysis;
  - `j2`, last scale for analysis;
  - `Jref`, reference scale under which eigenvalues have same bias;
  - `wtype`, kind of weighting used in the linear regressions:
    - 0 for no weighting  (i.e., uniform weights);
    - 1 for 1/nj weights with nj the wavelet sample size at scale j (suitable for fully Gaussian data);
    - 2 to use variance of the estimates;
  - `NB`, number of bootstrap resampling;
  - `LB`, number of blocks for the bootstrap resampling.

The main parameters contained in the structures `est` and `estbc` returned by OFBM_estimBC_BS are:
  - `est.hU`, matrix of univariate-like self-similarity exponent and cross-exponent estimates;
  - `est.h`, classical multivariate self-similarity exponent estimates;
  - `estbc.h`, bias corrected multivariate self-similarity exponent estimates.
    
The clustering of the scaling exponents need to run the estimation with adapted parameters:
```
% parameters of the estimation for clustering (bootstrap estimates are needed for the pairwise tests)
paramsEst.NB = 500; paramsEst.LB = 2*params.Nwt; 
% return self-similarity exponent estimation
[est,estbc] = OFBM_estimBC_BS(data,paramsEst) ;

% testing procedure
alpha = 0.05;
estT = OFBM_estimBC_BS_test(estbc,alpha);
% cluster the self-similarity exponents
[nbcluster,cluster] = successiveTestClustering(estT.decsortHocpw);
```
The parameters `alpha` of OFBM_estimBC_BS_test is the significance level of the multiple hypothesis test.


Another clustering method is available:
```
% cluster the self-similarity exponents
[nbcluster,cluster] = successiveTestClustering(estT.decsortHocpw_v2);
```

Examples with simulated ofBm can also be found in the `examples` folder.
