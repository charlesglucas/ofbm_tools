function [rout,K,funch] = MapAutoCovar(rin,name,params,forward)
% Mapping between nonGaussian and Gaussian auto-covariance sequences, rY
% and rX, respectively
%
% The non-Gaussian series is given by
% 1) Y[n]=f(X[n]), where X is a Gaussian series with covariance rX and mean
% zero
% or
% 2) Y[n]=f(X1[n],...,XK[n]), where X1,...,XK are iid Gaussian series each
% with covariance rX
%
% Usage:
%   [rout,K,funch] = MapAutoCovar(rin,name,params,forward)
% Inputs:
%   rin       the covariance sequence to transform, see parameter forward
%   "name" and "params" are used to choose the transform f, 
%           name is a "string" and "params" is a vector:
%             name    | params | description
%            'Chi2'    [d]       Y[n] = sum_{k=1}^d X_k[n]^2
%            'Exp'     [mu]      Y[n] = mu/2 (X_1[n]^2 + X_2[n]^2)
%            'Erlang'  [a,b]     Y[n] = b/2 sum_{k=1}^{2a} X_k[n]^2
%            'LN'      [mu,s2]   Y[n] = exp(sqrt(s2)X[n] + mu)
%            'Pareto'  [a,b]     Y[n] = b*exp((X_1[n]^2 + X_2[n]^2)/(2a))
%            'Uniform' [a,b]     Y[n] = a+(b-a)*exp(-0.5(X_1[n]^2 + X_2[n]^2))
%            'Laplace' [a]       Y[n] = a/2 (X_1[n]^2 - X_2[n]^2 + X_3[n]^2 - X_4[n]^2)
%            'ALaplace' [ap,am]       Y[n] = ap/2 (X_1[n]^2 + X_2[n]^2) -  am/2(X_3[n]^2  + X_4[n]^2)
%            'Gauss    [mu,s2]   Y[n] = sqrt(s2)*X[n]
%   forward  flag to choose between the mappings
%               i) backward, rY->rX, in this case rin=rY and rout=rX
%               ii) forward, rX->rY, in this case rin=rX and rout=rY
%            set to 1 or true to choose forward mapping
%            if parameter is omitted, the function will return backward
%            mapping
%
% Output
%   rout    the mapping of rin, see parameter forward
%   K       dimensionality of the Gaussian series X
%   funch   function handle for the mapping f
%
% Conditions on rX and parameters:
%             'Chi2' -- rX(k)>=0, d is a postitive integer
%             'Exp'  -- rX(k)>=0, mu>0
%             'LN'   -- rX(k)>-mu^2, s2>=0
%PEND: Finish conditions...
%
% Copyright 2009, Hannes Helgason, ENS de Lyon

if nargin<4,
    % default is backward mapping, i.e., from rY to rX
    forward = false;
end

switch upper(name),
    case 'CHI2'
        d = params(1);
        if sum(rin<0)>0,
            error('All entries in input covariance sequence need to be non-negative!');
        end
        if d<=0,
            error('Parameter "d=params(1)" is not positive!');
        end
        if forward,
            % mapping is rX->rY
            rout = 2*d*rin.^2;
        else
            % mapping is rY->rX
            rout = sqrt(rin/(2*d));            
        end
        K =  d; % dimensionality of Gaussian series
        funch = @(xx) sum(xx.^2); % function handle for the mapping X->Y
    case 'EXP'
        mu = params(1);
        if sum(rin<0)>0,
            error('All entries in input covariance sequence need to be non-negative!');
        end
        if mu<=0,
            error('Parameter "mu=params(1)" is not positive!');
        end
        if forward,
            % mapping is rX->rY
            rout = mu^2*rin.^2;
        else
            % mapping is rY->rX
            rout = sqrt(rin)/mu;
        end        
        K =  2; % dimensionality of Gaussian series
        funch = @(xx) mu/2*sum(xx.^2); % function handle for the mapping X->Y
    case 'ERLANG'
        a = params(1);
        b = params(2);
        if forward,
            % mapping is rX->rY
            rout = (a*b^2)*rin.^2;
        else
            % mapping is rY->rX
            rout = sqrt(rin/(a*b^2));
        end        
        K = 2*a; % dimensionality of Gaussian series
        funch = @(xx) b/2*sum(xx.^2); % function handle for the mapping X->Y
    case 'LN'
        mu = params(1);
        s2 = params(2);
        tmp = exp(-s2-2*mu);
        if forward,
            % mapping is rX->rY
            rout = (exp(s2*rin)-1)/tmp;
        else
            % mapping is rY->rX
            rout = log(1+tmp*rin)/s2;
        end
        K = 1; % dimensionality of Gaussian series
        funch = @(xx) exp(sqrt(s2)*xx + mu); % function handle for the mapping X->Y
    case 'GAUSSIAN'
        mu = params(1);
        s2 = params(2);
        if forward,
            % mapping is rX->rY
            rout = s2*rin;
        else
            % mapping is rY->rX
            rout = rin/s2;
        end
        K = 1; % dimensionality of Gaussian series
        funch = @(xx) sqrt(s2)*xx + mu; % function handle for the mapping X->Y  
	case 'PARETO'
        a = params(1);
        b = params(2);
        if forward,
            % mapping is rX->rY
            rout = a^2*b^2/(a-1)^2 * (rin.^2)./((a-1)^2 - rin.^2);
        else
            % mapping is rY->rX
            tmp = (a-1)^2/(a^2*b^2) * rin;
            rout = sqrt((a-1)^2 * tmp./(tmp + 1));
        end
        K =  2; % dimensionality of Gaussian series
        funch = @(xx) b*exp(sum(xx.^2)/(2*a));% function handle for the mapping X->Y         
	case 'UNIFORM'
        a = params(1);
        b = params(2);
        if forward,
            % mapping is rX->rY
            rout = (b-a)^2 * rin.^2./(16 - 4*rin.^2);
        else
            % mapping is rY->rX
            tmp = 4/(b-a)^2 * rin;
            rout = sqrt(4 * tmp./(tmp + 1));
            %rout = 4*sqrt( r./((b-a)^2+4*r) );
        end
        K = 2; % dimensionality of Gaussian series
        funch = @(xx) a+(b-a)*exp(-sum(xx.^2)/2);% function handle for the mapping X->Y         
	case 'LAPLACE'
        a = params(1);
        if forward,
            % mapping is rX->rY
            rout = 2*a^2*(rin).^2;
        else
            % mapping is rY->rX
            rout = sqrt(rin/(2*a^2));
        end
        K = 4; % dimensionality of Gaussian series
        % Check whether a/4 ot a/2 below
        funch = @(xx) a/2*(xx(1,:).^2 - xx(2,:).^2 + xx(3,:).^2 - xx(4,:).^2); % function handle for the mapping X->Y     
     case 'ALAPLACE'
        ap = params(1); am = params(2); a2 = (ap^2 + am^2)/2 ;
        if forward,
            % mapping is rX->rY
            rout = 2*a2*(rin).^2;
        else
            % mapping is rY->rX
            rout = sqrt(rin/(2*a2));
        end
        K = 4; % dimensionality of Gaussian series
        funch = @(xx) ap/2*(xx(1,:).^2 + xx(2,:).^2) - am/2*(xx(3,:).^2 + xx(4,:).^2) ; % function handle for the mapping X->Y         
    otherwise
        error(sprintf('Unsupported mapping: %s',name));
end
