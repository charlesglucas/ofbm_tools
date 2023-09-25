function rcirc = embedCrossCov(rxy,ryx)
% Circularization of cross-covariance
%
% Usage:
%   rcirc = embedCrossCov(rxy,ryx)
% Inputs:
%   The cross-covariance to circularize,
%     rxy    corresponds to r_XY[k]=EX[k]Y[0], k>=0
%     ryx    corresponds to r_YX[k]=EY[k]X[0]=r_XY[-k], k>=0
% Output
%   rcirc  the circularized cross-covariance
%
% Copyright (c) Hannes Helgason, ENS Lyon, 2009

N = length(ryx);
tmp = ryx(2:N);
rcirc = [rxy(1:(N-1)) fliplr(tmp)];
