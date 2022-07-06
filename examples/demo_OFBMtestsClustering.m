%% estimate and cluster multivariate self-similarity exponents
% CGL, October 2021.

clc
clear all
close all
format compact

addpath('../include/')
addpath(genpath('../WLBMF_tool/'))
load('../data/result_estimbc_sizeH6.mat')

%% Estimation and Bootstrap
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

%% Testing procedure
alpha = 0.05;
%paramsTest.P = size(data,1); paramsTest.NB = paramsEst.NB;
estT = OFBM_estimBC_BS_test(estbc,alpha);

disp(['H =                        [',sprintf(' %.1f ',params.H),']'])
disp(['H estimate =               [',sprintf(' %.2f ',estbc.h),']'])

%% Clustering based on sorted pairwise tests
[nbcluster, cluster] = successiveTestClustering(estT.decsortHocpw);
disp(['Sorted tests:              clusters = [',num2str(cluster),'], ',num2str(nbcluster),' clusters'])

%% Modified self-similarity exponent values
disp(['Adapted scaling exponents: [',sprintf(' %.2f ',averagedClusters(estT.h,cluster)),']'])

%% Clustering based on alternative sorted pairwise tests
[nbcluster_v2, cluster_v2] = successiveTestClustering(estT.decsortHocpw_v2);
disp(['Alternative sorted tests:  clusters = [',num2str(cluster_v2),'], ',num2str(nbcluster),' clusters'])

%% Modified self-similarity exponent values
disp(['Adapted scaling exponents: [',sprintf(' %.2f ',averagedClusters(estT.h,cluster_v2)),']'])