function U = BlockDiagMultiplication(B,W)
% Function for multiplying block-diagonal matrix diag([B1 B2 ... BM]),
% (Bk of size P-by-P) a vector [W1^T W2^T ... Wm^T]^T
% (Wk of column vector of length P). That is, calculates Bk*Wk
% 
% Usage:
%   U = BlockDiagMultiplication(B,W)
% Inputs:
%   B  matrix of size P-by-PM where
%      B(:,(k-1)*P+(1:P)) stores block Bk (of size P-by-P)
%   W  matrix of size P-by-M where
%      W(:,k) stores the column vector Wk
% Output:
%   U  matrix of size P-by-M where U = [U1 U2 ... UM], with Uk=Bk*Wk
%
% Copyright (c) 2010, Hannes Helgason, ENS de Lyon

%   Note: This M-file is intended to document the algorithm.
%   If a MEX file for a particular architecture exists,
%   it will be executed instead, but its functionality is the same.

% Implementing this in mex is orders of magnitude faster than using a for
% loop in Matlab (in the case M>>P, for example M=100,000, P=2)

[P,M] = size(W);

U = zeros(P,M); % using a matrix structure to store the vectors U1, U2,...,UP in rows
for m=1:M,
	Bm = B(:,(1:P)+(m-1)*P);
	
tmp = Bm*W(:,m);
	U(:,m) = tmp;
end