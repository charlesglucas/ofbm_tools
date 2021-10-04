# OFBM tools
===

## Description
SUGAR DMS is an image recontruction and contour detection matlab package. It minimizes the Discrete Mumford-Shah (D-MS) functional, which enforces constraints related to smoothing over image and sparsity over contours to obtain a piecewise smooth image from observed noisy images, with an automatic selection of the regularization parameters based on a Stien-like strategy.

## Requirements
The following matlab package is required: [GRANS0](https://gitlab.com/timmitchell/GRANSO/).

## References
  - [Lucas et al., 2021] - Preprint
  - [Foare et al., 2019](https://hal.archives-ouvertes.fr/hal-01782346/document)
  - [Deledalle et al., 2014](https://arxiv.org/pdf/1405.1164)
  
## Quick start
The basic syntax to run SUGAR D-MS is as follows:

```
% return optimal hyperparameters of D-MS
[Lambda,~] = bfgs_sugar_dms(image);

% return D-MS image reconstruction u and contour e estimates
[u,e,~] = DMS_2D(image,Lambda(1),Lambda(2));
```


The main parameters to take into account are:

  - `R`, the number of realizations of the Monte Carlo vector;
  - `sigma`, the noise level which is estimated by default.
    
Here is an example with non-default parameters:
```
param.R = 5; param.sigma = 0.1;
[Lambda,~] = bfgs_sugar_dms(image, param);
[u,e,~] = DMS_2D(image,Lambda(1),Lambda(2));
```
An example with simulated images can also be found in the `example` folder.
