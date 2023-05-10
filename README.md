# <div align="center">OFBM Tools </div> 

## References

> **[Charles-Gérard Lucas](https://perso.ens-lyon.fr/charles.lucas), [Herwig Wendt](https://www.irit.fr/~Herwig.Wendt/), [Patrice Abry](https://perso.ens-lyon.fr/patrice.abry), [Gustavo Didier](http://www2.tulane.edu/~gdidier/),**
*Multivariate time-scale bootstrap for testing the equality of selfsimilarity parameters,* 
Colloque Francophone de Traitement du Signal et des Images (GRETSI). [Download](https://hal.archives-ouvertes.fr/hal-03735529)

> **[Charles-Gérard Lucas](https://perso.ens-lyon.fr/charles.lucas), [Patrice Abry](https://perso.ens-lyon.fr/patrice.abry), [Herwig Wendt](https://www.irit.fr/~Herwig.Wendt/), [Gustavo Didier](http://www2.tulane.edu/~gdidier/),**
*Counting the number of different scaling exponents in multivariate scale-free dynamics: Clustering by bootstrap in the wavelet,* International Conference on Acoustics, Speech, & Signal Processing (ICASSP). [Download](https://hal.archives-ouvertes.fr/hal-03735481/document)

## Description
OFBM Tools is a multivariate self-similarity analysis matlab package. It permits to estimate and count the values of scaling exponents $H_1,\ldots,H_M$ across $M$-variate time series modeled by $M$-variate fractional Brownian motion ($M$-fBm) $\underline{B}_{\underline{\underline{H}},\Sigma}$. 

![alt text](http://perso.ens-lyon.fr/charles.lucas/images/multivariateHestim.svg)

The estimation is based on linear regressions of wavelet eigenvalues across scales. A de-biased estimator based on the computation of wavelet spectra from equal number of wavelet coefficients is designed. To count the scaling exponents, a clustering strategy is derived from parametric pairwise hypothesis tests between successive ordered estimates: two bootstrap-based estimation methods are proposed for the computation of the test parameters.

## Recommendation
This toolbox is designed to work with [**Matlab 2020b**](https://fr.mathworks.com/products/new_products/release2020b.html).

## Quick start

<details open>
  <summary><strong>Estimation</strong></summary>

The basic syntax to run `OFBM_estimBC_BS` is as follows:

```
paramsEst.Nwt = 2 ; paramsEst.FigNum = 1 ; paramsEst.wtype = 1 ;
paramsEst.j1 = 8; paramsEst.j2 = 11; paramsEst.Jref = paramsEst.j2 ; 
paramsEst.NB = 0; paramsEst.LB = 0;
[est,estbc] = OFBM_estimBC_BS(data,paramsEst);
```
The input `data` is a $M \times N$ matrix with $M$ the number of components and $N$ the sample size.
  
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
    - 2 to use the variance of the estimates;
  - `NB`, number of bootstrap resampling;
  - `LB`, number of blocks for the bootstrap resampling.

The main parameters contained in the structures `est` and `estbc` returned by OFBM_estimBC_BS are:
  - `est.hU`, matrix of univariate-like self-similarity exponent and cross-exponent estimates;
  - `est.h`, classical multivariate self-similarity exponent estimates;
  - `estbc.h`, bias corrected multivariate self-similarity exponent estimates.
</details>

<details open>
  <summary><strong>Clustering</strong></summary>

The count of the self-similarity exponents needs to run `OFBM_estimBC_BS` with adapted parameters `paramsEst` for the bootstrap procedure:
```
paramsEst.NB = 500; paramsEst.LB = 2*params.Nwt; 
[est,estbc] = OFBM_estimBC_BS(data,paramsEst) ;
```

The routine `OFBM_estimBC_BS_test` gives decisions $d_{\alpha}^{(m)}$ for hypothesis $H_m=H_{m+1}$ with a false discovery rate $\alpha$. These decisions naturally separate the estimates in different clusters.
<p align="center">
  <img width="500" src="http://perso.ens-lyon.fr/charles.lucas/images/naiveClustering.svg" style="max-width: 100%;">
  <img align="center" width="1500" src="https://github.com/HoangTrieuVy/Deep-Unfolding-Networks-for-Image-Denoising/raw/main/examples/results_unfolded_ScCP_10081_noisy.jpg" style="max-width: 100%;">
</p>

The clustering strategy, with a false discovery rate `alpha` for the multiple hypothesis test and the bootstrap-based test parameter estimation described in [Lucas et al., ICASSP 2022](https://hal.archives-ouvertes.fr/hal-03735481/document), can be run as follows:
```
alpha = 0.05; estT = OFBM_estimBC_BS_test(estbc,alpha);
[nbcluster,cluster] = successiveTestClustering(estT.decsortHocpw);
```

Another bootstrap-based test parameter estimation, described in [Lucas et al., GRETSI 2022](https://hal.archives-ouvertes.fr/hal-03735529), is available for the clustering method:
```
[nbcluster,cluster] = successiveTestClustering(estT.decsortHocpw_v2);
```

Examples with synthetic ofBm can also be found in the `examples` folder.
</details>
