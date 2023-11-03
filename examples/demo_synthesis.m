%% synthetize OFBM 
% CGL, September 2023.

clc
clear all
close all
format compact

addpath('../includeSynthesis/')

%% synthesis parameters
params.H = [.6 .6 .6 .8 .8 .8];
P = length(params.H) ; params.P = P ; 
%params.W = RandOrthMat(P);
params.W = RandomMixing(P,0.5); params.W = MixingNorm(params.W) ;
rho = .5 ; params.Correl = toeplitz([1 rho*ones(1,P-1)]);  params.CorrDiag = 0 ; 
params.nbsamples = 2^16 ;
params.FBM = 1;
params.Sigma2 = ones(1,length(params.H)); 
params.Nwt = 2 ;

%% generate data

if params.CorrDiag == 0  
    params.nbcopies = 0;
    [~,~,Rxapprox,approxFlag,BmArray,RX,K,funch] = OFBMCovSynthesis(params) ;  
else
    K=0;RX=0;BmArray=0;funch='' ; 
end
[data,dataNoMix] = OFBMSynthesis(K,RX,BmArray,funch,params) ; 
