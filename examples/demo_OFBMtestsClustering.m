%% estimate and cluster multivariate self-similarity exponents
% CGL, October 2021.

clc
clear all
close all
format compact

load('data/result_estimbc_sizeH6.mat')

%% Estimation and Bootstrap
params.Nwt = 2 ;
params.j1 = 6 ;
params.j2 = 10 ;
params.FigNum = 10 ;
params.wtype = 1 ;
params.Jref = params.j2;
params.NB=500;
params.LB=params.Nwt;
params.Bcorr=0;

[est,estbc] = OFBM_estimBC_BS(data,params) ;

%% Testing procedure
alpha = 0.05;
estT = OFBM_estimBC_BS_test(estbc,alpha,params);

%% Clustering based on sorted pairwise tests
[nbcluster, cluster] = successiveTestClustering(estT.decsortHocpw);
disp(['Sorted tests: clusters = [',num2str(cluster),'], ',num2str(nbcluster),' clusters'])