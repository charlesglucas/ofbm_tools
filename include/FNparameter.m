function [mu,sigma] = FNparameter(y,xinit)
% returns Folded Normal Parmeters from samples y
%
% Charles-Gérard Lucas, ENS Lyon, 2021

if nargin<2, xinit = 0.9*sqrt(mean(y.^2)); end
N = length(y);
fun = @(x) f(x,y,N);
x = fzero(fun,xinit);
mu = abs(x);
sigma = sqrt(max(0,mean(y.^2)-mu^2));
end

function z = f(x,y,n)
sigma2 = mean(y.^2) - x.^2;
z = sum(y.*(1-exp(2*x*y/sigma2))./(1+exp(2*x*y/sigma2))) + n*x;
end
