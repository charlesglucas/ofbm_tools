function rcirc = embedAutoCov(r)
% Circularization of auto-covariance
%
% Usage:
%   rcirc = embedAutoCov(r)
% Inputs
%   r      the auto-covariance sequence to circularize
% Output
%   rcirc  the circularized auto-covariance
%
% Copyright (c) Hannes Helgason, ENS Lyon, 2009

N = length(r);
tmp = r(2:N-1);
rcirc = [r(1) tmp r(N) fliplr(tmp)];
