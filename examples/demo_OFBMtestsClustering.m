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
paramsEst.Nwt = 2 ;
paramsEst.j1 = 8 ;
paramsEst.j2 = 11 ;
paramsEst.FigNum = 1 ;
paramsEst.wtype = 1 ;
paramsEst.Jref = paramsEst.j2;
paramsEst.NB = 500;
paramsEst.LB = 2*paramsEst.Nwt;

[est,estbc] = OFBM_estimBC_BS(data,paramsEst) ;

disp(['H =                        [',sprintf(' %.1f ',params.H),']'])
disp(['H estimates =               [',sprintf(' %.2f ',estbc.h),']'])

%% Testing if all exponents are equal or not
alpha = 0.05;
testChi2 = BSChi2test(estbc,alpha);
disp(['P-value = ',sprintf('%.2f',testChi2.pval)])
disp(['Decision = ',num2str(testChi2.dec)])

%% Clustering based on M-1 half normal pairwise tests
alpha = 0.05;
testHN = BSHalfNormalTest(estbc,alpha);
[nbclusterHN,clusterHN] = successiveTestClustering(testHN.decsortHocpw);
disp(['Sorted tests:              clusters = [',num2str(clusterHN),'], ',num2str(nbclusterHN),' clusters'])
disp(['Adapted scaling exponents: [',sprintf(' %.2f ',averagedClusters(est.h,clusterHN)),']'])

%% Clustering based on M-1 folded normal pairwise tests
alpha = 0.05;
testFN = BSFoldedNormalTest(estbc,alpha);
[nbclusterFN, clusterFN] = successiveTestClustering(testFN.decsortHocpw);
disp(['Alternative sorted tests:  clusters = [',num2str(clusterFN),'], ',num2str(nbclusterFN),' clusters'])
disp(['Adapted scaling exponents: [',sprintf(' %.2f ',averagedClusters(est.h,clusterFN)),']'])