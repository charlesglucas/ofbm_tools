function [Rapprox,approxFlag,BmArray] = InitStepMultivarGaussBestApprox(R,N,genTR)
% Initialize structures for the generationt multivariate stationary Gaussian 
% processes provided the spectrum.
% Uses the alternative construction (\tilde{X}=F^* U) introduced
% in the multivariate non-Gaussian synthesis paper.
%
% If R is not a valid covariance structure (satisfying
% nonnegative definiteness) it does an approximation (see non-Gaussian
% paper)
%
% Usage:
% [Rapprox,approxFlag,BmArray] = InitStepMultivarGaussBestApprox(R,N,genTR)
%
% Inputs:
%   R, N, genTR  see documentation given in GenMultivarGaussBestApprox
% Outputs:
%   X and Xiid are two independent realizations of the process, stored
%   as matrices of size P-by-N. X(k,:) is a realization of Xk, i.e., 
%   component k, and is a vector of length N. Same holds for Xiid.
%
%   approxFlag  true if approximation was needed, otherwise false
%
%   BmArray     this structure can be passed to
%               SynthStepMultivarGaussBestApprox to generate samples of
%               the process, e.g., for Monte Carlo simulations
% 
% Copyright (c) Hannes Helgason, ENS Lyon, 2010

P = size(R,1); % number of coordinates

if P>1,
    % Perform checks for multivariate processes
    if nargin==3,
        if ~genTR & isempty(R{2,1}),
            warning('genTR is set to false but R{2,1} is empty. Will assume symmetric cross-covariance.'); 
        end
    elseif nargin<3,
        genTR = false;
    end

    if genTR | isempty(R{2,1}),
        % Using symmetric cross-covariance (time reversibility)
        asymm = false;
    else
        asymm = true;
    end
else
    % the process is monovariate
    asymm = false;
end
    
% circular embedding of covariance
Rcirc = circembed(R,asymm);

% Calculate auto- and cross-spectrum
Lambda = cell(P,P); % initialize cell structure
for k=1:P,
  for j=k:P,
    Lambda{k,j} = fft(Rcirc{k,j});
  end
end

% and the length of the circularly embedded covariance sequences
M = length(Lambda{1,1});

% For returning the resulting covariance sequence
LambdaApprox = cell(P,P);
for k=1:P,
    for j=1:P,
        LambdaApprox{k,j} = zeros(1,M);
    end
end

% initialize
Gm = zeros(P,P);
Gupper = zeros(P,P);
Bmarray = zeros(P,M*P);
approxFlag = false; % flag which tells whether approximation was needed

for m=1:M,
  % fill in upper-half of the matrix Gm
  for k=1:P,
    for j=(k+1):P,
%      Gupper(k,j) = 2*Lambda{k,j}(m);
      Gupper(k,j) = Lambda{k,j}(m);
      %PEND: Multiplication by two is done in circembed, should be better
      %to do here
    end
  end
  % the lower-half of Gm is the conjugate transpose of its upper-half
  Gm = Gupper + Gupper';
  
  % fill in the diagonal entries
  for k=1:P,
    % The auto-spectrum should be real-valued, taking real-part in case of
    % numerical errors
    tmp = real(Lambda{k,k}(m)); % taking real part in case of numerical errors
    Gm(k,k) = 2*tmp;
%    Gm(k,k) = tmp;
  end

  % At this point, Gm could have negative eigenvalues.
  % Perform Schur decomposition of this Hermitian matrix
  % and construct the nonnegative definite matrix which best approximates
  % Gm under the Frobenius norm if needed
  [Q,S] = schur(Gm); % what about SVD?
  diagS = real(diag(S));

  diagStilde = zeros(P,1);
  for k=1:length(diagS),
	if diagS(k) >= 0,
        diagStilde(k) = diagS(k);
    else
      % approximation is needed, replace negative value with zero
	  approxFlag = true;
      diagStilde(k) = 0;
    end
  end

  % Now we could construct the approximation Gmtilde = Q*Stilde*(Q'), 
  % but instead we directly construct Bm in Gmtilde=Bm*Bm'
  tmp = sqrt(diagStilde);
  sqrtStilde = diag(tmp);
  Bm = Q*sqrtStilde;
  BmArray(:,(1:P)+(m-1)*P) = Bm;
  
  % At this step we could generate Um this way:
  %  W = (randn(P,1) + i*randn(P,1))/sqrt(2);
  %  tmp = Bm*W;
  %  U(:,m) = tmp;
  % This is now done in SynthStepMultivarGaussBestApprox
  
  % Calculate Gmtilde to be able to get the resulting covariance sequence
  Gmtilde = Bm*(Bm');
  for k=1:P,
    for j=1:P,
        LambdaApprox{k,j}(m) = Gmtilde(k,j)/2;
    end
  end
end

% calculate the resulting covariance structure
Rapprox = cell(P,P);
nR = length(R{1,1});
for k=1:P,
	for j=k:P,
        tmp = ifft(LambdaApprox{k,j});
        Rapprox{k,j} = tmp(1:N);
        tmpFlip = fliplr(tmp((N+1):end));
        % fill in for R(k,j) with j<k (cross-covariance does not need to be
        % symmetric)

        Rapprox{j,k} = [tmp(1) tmpFlip(1:end)];
      
    end
end