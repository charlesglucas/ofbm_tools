%% try conquere-and-divide bias correction
% CGL, October 2021.

clc
clear all
close all 
format compact

addpath(genpath('./'))
load('data/result_estimbc_sizeH4.mat')

%% Estimation
paramsEst = params;
paramsEstEst.Nwt = 2 ;
paramsEst.j1 = 5;
paramsEst.j2 = 10;
paramsEst.Jref = paramsEst.j2 ;
paramsEst.FigNum = 10 ;
paramsEst.wtype = 1 ;
paramsEst.NB = 0;
paramsEst.LB = 0;

[est,estbc] = OFBM_estimBC_BS(data,paramsEst) ;
disp(['Classical H estimate:        [',num2str(est.h),']'])
disp(['Bias corrected H estimate:   [',num2str(estbc.h),']'])