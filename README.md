OFBM tools
===

## Description
OFBM tools is a multivairate sel-similarity analysis matlab package. It permits to estimate and group the values of scaling exponents across multivariate time series modeled by an operator fractional Brownian motion (ofBm).

## Requirements
The following matlab package is required: [WLBMF](https://www.irit.fr/~Herwig.Wendt/software.html).

## References
  - [Lucas et al., 2021](https://www.irit.fr/~Herwig.Wendt/data/LucasEUSIPCO2021.pdf)
  
## Quick start
The basic syntax to run OFBM tools is as follows:

```
paramsEst = params;
paramsEst.Nwt = 2 ;
paramsEst.j1 = 5;
paramsEst.j2 = 10;
paramsEst.Jref = paramsEst.j2 ;
paramsEst.FigNum = 10 ;
paramsEst.wtype = 1 ;
paramsEst.NB = 0;
paramsEst.LB = 0;
% return self-similarity exponent estimates
[est,estbc] = OFBM_estimBC_BS(data,paramsEst) ;
```



The main parameters to take into account in the structure `params` are:

  - `R`, the number of realizations of the Monte Carlo vector;
  - `Jref`: reference scale under which several wavelet spectra are computed with the same number wavelet coeficients
  - `j2`: last scale for analysis
  - `j1`: first scale for analysis
  - `wtype`, the kind of weighting used in the linear
                   regressions:
                       0 -  no weigthing  (ie uniform weights)
                       1 -  1/nj weights  (suitable for fully Gaussian data)
                       2 -  use variance of the estimates
  - `NB`: number of bootstrap resampling
  - `LB`: number of blocks for the bootstrap resampling
  

    
Here is an example with non-default parameters:
```
param.R = 5; param.sigma = 0.1;
[Lambda,~] = bfgs_sugar_dms(image, param);
[u,e,~] = DMS_2D(image,Lambda(1),Lambda(2));
```
An example with simulated images can also be found in the `example` folder.
