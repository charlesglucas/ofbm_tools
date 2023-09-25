function R = CovarMfGn_sigma_v2(N,h,sigmae)
% Multivariate fractional Gaussian noise covariance generator
%
% Usage:
%   R = CovarMfGn(N,h,W,A)
% Inputs:
%   N     length of the sequence
%   h     P-by-1 vector of the diagonal entries in the matrix D below
%         NOTE: at the moment we consider the entries to be real and
%         non-zero
%   sigmae  covariance matrix at instant 1
%
% Output
%   R     generated covariance structure, a cell structure of size P-by-P.
%
% Comments:
%   -- W is an invertible matrix, assumed diagonal see below
%   --Operator/Multivariate Fractional Brownian Motion (OFBM/MFBM)--
%   Consider a multivariate process B_H(t) = (B1(t),...,Bp(t))^T,
%   satisfying the self-similarity relation B_H(ct) = c^H B_H(t), c>0,
%   where H is a matrix. The equality should be understood in the sense
%   that the rhs and lhs have the same distribution. 
%   Recall the definition of matrix exponents, so that for a real number
%   c>0, we have
%     c^H=exp(H log(c)) = sum_{k=0}^\infty H^k (log(c))^k/k!
%   If B_H(t) is time-reversible (i.e., B_H(t) and B_H(-t) have the same
%   distribution) then
%     E B_H(t)B_H(s)* = 
%         1/2*{ -|t-s|^H G |t-s|^H* + |t|^H G |t|^H* + |s|^H G |s|^H* }
%     where G=E B_H(1)B_H(1)*
%
%   --Multivariate Fractional Gaussian Noise (MfGn)--
%   Let X_n = B_H(n+1) - B_H(n)
%   Then we can show
%
%     EX_n X_0* = 
%        {|n+1|^H G |n+1|^H* + |n-1|^H G |n-1|^H* - 2 |n|^H G |n|^H* }
%
%   Note that in the univariate case, G and H are real numbers, and
%     EX_n X_0 = 
%        G {|n+1|^(2H) + |n-1|^(2H) - 2 |n|^(2H) }
%
%   Assume we can diagonalize H, s.t., H=W*D*inv(W), we get
%     a^H = W a^D inv(W)
%   so that
%     EX_n X_0* = 
%        W {|n+1|^D inv(W) G inv(W)* |n+1|^D* + 
%           |n-1|^D inv(W) G inv(W)* |n-1|^D* +
%           |n|^D inv(W) G inv(W)* |n|^D*} W*
%
%   Write F=inv(W) G inv(W)*, and F_n = |n|^D F |n|^D*
%   Then the change of basis Y_n = inv(W)X_n gives
%     EY_n Y_0* = F_{n+1} + F_{n-1} - 2 F_n
%   
%   The matrix A comes from the spectral representation (see Didier and
%   Pipiras (2009)):
%     B_H(t) = 
%     \int_R (e^{itx}-1)/(ix) (x_+^{-(H-1/2)}A + x_-^{-(H-1/2)}bar{A})B(dx)
%   where x_+ = max{x,0}, x_- = max{-x,0}, bar{A} is the complex conjugate
%   of A and B(dx) is a suitable complex-valued Brownian motion.
%
% Copyright 2009, Hannes Helgason, ENS Lyon
% Modified on 05-20-2015, JF



P = length(h); % dimensionality



% Calculate G=E B_H(1)B_H(1)*
G = sigmae;

% Variable storing covariance, RR(:,:,n+1)=EX_n X_0*
RR = zeros(P,P,N);

% Calculate covariance for n=0 
F = G;
Fnmin1 = F; % since |-1|^D = Id
Fn = 0; % since |0|^D = 0
Fnplus1 = F; % since |1|^D = Id

% Note: EX_nX_0*= W(EY_nY_0*)W*, where EY_nY_0*=Bnplus1 - 2*Bn + Bnmin1
RR(:,:,1) = (Fnplus1 - 2*Fn + Fnmin1);

% Calculate covariance for n>0 
for n=1:(N-1),
  Fnmin1 = Fn;
  Fn = Fnplus1;
  tmp = (n+1).^h;
  Lambda = (tmp.') * tmp;
  Fnplus1 = Lambda.*F;
  RR(:,:,n+1) = (Fnplus1 - 2*Fn + Fnmin1);
end

% write the result in a cell array of size P-by-P

R = cell(P,P);
for k=1:P,
    for j=k:P,
        % Due to time-reversibility we can leave R{k,j} empty for k>j
        R{k,j} = (1/2)*squeeze(RR(k,j,:)).'; % take transpose to get a row vector %/!\ Rajout de 1/2
    end
end




