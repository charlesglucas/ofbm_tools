function [test] = BSChi2test(est,alpha)

% returns the chi-squared testing procedure 
%
%  Input:     est:     structure of estimation method 
%             alpha:   significance level of the tests
%
%  Output:    test:    structure of folded normal testing procedures
%                       - T : chi2 statistic
%                       - dec : decision
%                       - pval : p-value
%
% Charles-Gérard Lucas, ENS Lyon, 2021

P = length(est.lambda);

% compute approximate chi2 test
test.invSigma = inv(cov(est.lambdaBS));
cL = est.lambda - mean(est.lambda);
test.T=cL*test.invSigma*cL';
test.dec = test.T >= chi2inv(1-alpha,P-1);
test.pval = 1 - chi2cdf(test.T,P-1);

end

