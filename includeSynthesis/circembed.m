function Rcirc = circembed(R,asymm)
% Circular embedding of covariance
%
% Usage:
%   Rcirc = circembed(R,asymm)
% Inputs
%   R      covariance structure to circularly embed
%   asymm   
% Output
%   Rcirc  the embedded covariance
%
% Copyright (c) Hannes Helgason, ENS Lyon, 2009

if nargin<2,
    asymm = false;
end

P = size(R,1);
Rcirc = cell(P,P);

% circularize auto-covariance sequences
for k=1:P,
    Rcirc{k,k} = embedAutoCov(R{k,k});
end

% circularize cross-covariance sequences
if asymm,
    % assuming asymmetric cross-covariance
    for k=1:P,
        for j=(k+1):P,
            % multiplying by two to get the correct cross-covariance
            % using complex-valued Gaussians in the simulation
            Rcirc{k,j} = 2*embedCrossCov(R{k,j},R{j,k});
        end
    end
else
    % assuming symmetric cross-covariance
    for k=1:P,
        for j=(k+1):P,
            % multiplying by two to get the correct cross-covariance
            % using complex-valued Gaussians in the simulation            
            Rcirc{k,j} = 2*embedCrossCov(R{k,j},R{k,j});
        end
    end
end