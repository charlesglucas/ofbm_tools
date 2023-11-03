%% cluster multivariate sorted self-similarity exponents 
% CGL, October 2021.

clc
clear all
close all
format compact

addpath('../include/')
load('../data/result_estimbc_sizeH6.mat')

%% Estimation and bootstrap
paramsEst.FBM = 1;
paramsEst.Nwt = 2;
paramsEst.j1 = 8;
paramsEst.j2 = 11;
paramsEst.FigNum = 1;
paramsEst.wtype = 1;
paramsEst.Jref = paramsEst.j2;
paramsEst.NB = 500;
paramsEst.LB = 2*paramsEst.Nwt;

[est,estbc] = OFBM_estimBC_BS(data,paramsEst);

disp(['H =                        [',sprintf(' %.1f ',params.H),']'])
disp(['H estimates =               [',sprintf(' %.2f ',estbc.h),']'])

%% Testing if all exponents are equal or not (Hartigan's test)
alpha = .05;
BSdipTest(estbc,alpha)
 
%% Couting modes (Silvermans's test)
alpha = .05;
[nbm,silv,pvalS] = SilvermanModeEstim(estbc.h,alpha);
clusters = kmeans(estbc.h',nbm);
disp([num2str(nbm),' clusters',': [',num2str(clusters'),']'])