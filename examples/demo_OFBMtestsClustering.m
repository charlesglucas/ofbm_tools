%% estimate and cluster multivariate self-similarity exponents
% CGL, October 2021.

clc
clear all
close all
format compact

load('data/result_estimbc_sizeH6.mat')

%% Estimation and Bootstrap
paramsEst.FBM = 1;
paramsEst.Nwt = 2 ;
paramsEst.j1 = 8 ;
paramsEst.j2 = 10 ;
paramsEst.FigNum = 0 ;
paramsEst.wtype = 1 ;
paramsEst.Jref = paramsEst.j2;
paramsEst.NB=500;
paramsEst.LB=2*paramsEst.Nwt;

[est,estbc] = OFBM_estimBC_BS(data,paramsEst) ;

%% Testing procedure
alpha = 0.05;
paramsTest.P = size(data,1); paramsTest.NB = paramsEst.NB;
estT = OFBM_estimBC_BS_test(estbc,alpha,paramsTest);

%% Clustering based on sorted pairwise tests
[nbcluster, cluster] = successiveTestClustering(estT.decsortHocpw);
disp(['Sorted tests: clusters = [',num2str(cluster),'], ',num2str(nbcluster),' clusters'])

%% Modify self-similarity exponent values
averagedClusters(estT.h,cluster)