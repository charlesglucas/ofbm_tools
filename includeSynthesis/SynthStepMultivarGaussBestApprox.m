function [X,Xiid] = SynthStepMultivarGaussBestApprox(BmArray,N)
% Generate multivariate stationary Gaussian process provided spectrum.
% Uses the alternative construction (\tilde{X}=F^* U) introduced
% in the multivariate non-Gaussian synthesis paper.
%
% Usage:
%   [X,Xiid] = SynthStepMultivarGaussBestApprox(BmArray,N)
%
% Inputs:
%   Bmarray, N  see documentation given in InitStepMultivarGaussBestApprox
% Outputs:
%   X and Xiid are two independent realizations of the process, stored
%   as matrices of size P-by-N. X(k,:) is a realization of Xk, i.e., 
%   component k, and is a vector of length N. Same holds for Xiid.
% 
% Copyright (c) Hannes Helgason, ENS Lyon, 2010

% find the number of components
P = size(BmArray,1);
M = size(BmArray,2)/P;

W = (randn(P,M) + i*randn(P,M))/sqrt(2);

U = BlockDiagMultiplication(BmArray,W);

% Generate each component
Xtilde = 1/sqrt(M)*fft(U,[],2); % taking DFT along rows
% Note: dividing by sqrt(M) to normalize correctly

% the first N row entries of the real and imaginary parts
% are independent realizations of the process (respectively)
X = real(Xtilde(:,1:N));
Xiid = imag(Xtilde(:,1:N));
