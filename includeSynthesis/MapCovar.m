function [Rout,K,funch] = MapCovar(Rin,nameList,paramList,isforward, varY)
% Mapping between nonGaussian and Gaussian covariance structures, RY
% and RX, respectively. RY
%
% The p-th coordinate of the non-Gaussian series is given by
% 1) Yp[n]=f(Xp[n]), where Xp is a Gaussian series with covariance RX with
% mean zero
% or
% 2) Yp[n]=f(X1p[n],...,XKp[n]), where X1,...,XK are iid multivariate 
% Gaussian series each with covariance RX
%
% Usage:
%   [Rout,K,funch] = MapCovar(Rin,name,params,forward)
% Inputs:
%   Rin       the covariance structure to transform, see parameter forward
%             stored in a P-by-P cell array where Rin{k,j} is a vector
%             corresponding to a covariance sequence r_{kj}[n], n=0,...,N
%   nameList  parameter for choosing the transforms (see also
%           documentation for MapAutoCovar.m). This parameter can be
%           either:
%           1) string of transform name, same tranform will be used for all
%              components
%           2) NOT IMPLEMENTED YET: string array of transform names, nameList{p} corresponds to
%              the name of the transform for component p
%               
%   paramList  parameters corresponding to the transforms in nameList. This
%           parameter can be either of the two:
%           1) vector, same parameters will be used for all components
%           2) NOT IMPLEMENTED YET: cell array of vectors, paramList{p} are the parameters used
%              for component p
%   isforward  flag to choose between the mappings
%               i) backward, RY->RX, in this case Rin=RY and Rout=RX
%               ii) forward, RX->RY, in this case Rin=RX and Rout=RY
%            set to 1 or true to choose forward mapping
%            if parameter is omitted, the function will return backward
%            mapping
%   
%
% Output
%   Rout    the mapping of Rin, see parameter forward
%           stored in a P-by-P cell array where Rin{k,j} is a vector
%           corresponding to a covariance sequence r_{kj}[n], n=0,...,N
%   K       number of independent multivariate Gaussian copies used in the 
%           transform of the Gaussian series X
%   funch   function handle for the mapping f
%
% Conditions on RX and parameters:
%PEND: Finish writing about conditions...
%
% Copyright 2010, Hannes Helgason, ENS de Lyon

if nargin<4,
    % default is backward mapping, i.e., from RY to RX
    isforward = false;
end

if ~isforward,
    % we are doing backward mapping, from rY->rX, so need to scale.
    Rin = NormalizeCovar(Rin,nameList,paramList);
end

if iscell(Rin),
    P = length(Rin);
else
    % we are in the univariate case
    P = 1;
    % put Rin into a cell array
    Rintmp = Rin;
    Rin = cell(1,1);
    Rin{1} = Rintmp;
end

Rout = cell(P,P);
for k=1:P,
    for j=1:P,
        if ~isempty(Rin{k,j}),
            [Rout{k,j},K,funch] = MapAutoCovar(Rin{k,j},nameList,paramList,isforward);

        end
    end
end



% normFactor = ones(1,P);    % normalization to get the variance for marginal,
%                             % if we are doing forward mapping,
%                             % the variance of the Gaussian marginal is 1
% if ~isforward,
%     % we are doing backward mapping, from rY->rX, so need to scale.
%     % Find the variance of unscaled marginal
%     for k=1:P,
%         if iscell(nameList),
%             error('nameList has to be a string: The transforms for each component have to be the same.');
%             %funchName = nameList{k};
%         else
%             funchName = nameList;
%         end
%         if iscell(paramList),
%             error('paramList has to be a string: The transforms for each component have to be the same.');            
%             funcParam = paramList{k};
%         else
%             funcParam = paramList;
%         end
%         varmarg = MapAutoCovar(1,funchName,funcParam,true);
%         normFactor(k) = sqrt(varmarg/Rin{k,k}(1));
%     end
% end
% 
% 
% Rout = cell(P,P);
% for k=1:P,
%     for j=1:P,
%         if ~isempty(Rin{k,j}),
%             % normalize the input covariance so it has the covariance of
%             % the chosen marginal
%             Rinkj_normalized = normFactor(k)*normFactor(j)*Rin{k,j};
%             [Rout{k,j},K,funch] = MapAutoCovar(Rinkj_normalized,nameList,paramList,isforward);
% %            [Rout{k,j},K,funch] = MapCrossCovar(Rinkj_normalized,funch1,param1,funch2,param2,forward);
%         end
%     end
% end
