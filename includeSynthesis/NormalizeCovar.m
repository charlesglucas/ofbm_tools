function Rout = NormalizeCovar(Rin,nameList,paramList)
% Normalizes input covariance structure so that the variance matches the one
% given in the transform defined by paramList
%
% Usage:
%   Rout = NormalizeCovar(Rin,nameList,paramList)
% Inputs:
%   Rin       the covariance structure to normalize
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
%
% Output
%   Rout    the normalized version of Rin
%           stored in a P-by-P cell array where Rout{k,j} is a vector
%           corresponding to a covariance sequence r_{kj}[n], n=0,...,N
%
% Copyright 2010, Hannes Helgason, ENS de Lyon

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

normFactor = ones(1,P);    % normalization to get the variance for marginal,
                            % if we are doing forward mapping,
                            % the variance of the Gaussian marginal is 1

% Find the variance of unscaled marginal
for k=1:P,
	if iscell(nameList),
        error('nameList has to be a string: The transforms for each component have to be the same.');
        %funchName = nameList{k};
    else
        funchName = nameList;
    end
    if iscell(paramList),
        error('paramList has to be a string: The transforms for each component have to be the same.');            
        funcParam = paramList{k};
    else
        funcParam = paramList;
    end
    varmarg = MapAutoCovar(1,funchName,funcParam,true);
    normFactor(k) = sqrt(varmarg/Rin{k,k}(1));
end


Rout = cell(P,P);
for k=1:P,
    for j=1:P,
        if ~isempty(Rin{k,j}),
            % normalize the input covariance so it has the covariance of
            % the chosen marginal
            Rout{k,j} = normFactor(k)*normFactor(j)*Rin{k,j};
        end
    end
end

if P == 1,
   % we are in the univariate case, return vector
   Rout = Rout{1,1};
end