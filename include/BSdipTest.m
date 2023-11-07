function [test] = BSdipTest(est,alpha)

% returns the Hartigan's dip test procedure
%
%  Input:     est:     structure of estimation method 
%             alpha:   significance level of the test  
%
%  Output:    test:    structure of half normal testing procedures
%                       - pval: p-value
%                       - dec: decision
%
% Charles-Gérard Lucas, ENS Lyon, 2023

NB = size(est.lambdaBS,1);
dip = HartigansDipTest(sort(est.h));
for r=1:NB
	dipBS(r) = HartigansDipTest(sort(est.lambdaBS(r,:)-mean(est.lambdaBS(r,:))));
end
ihi = round((1-alpha)*NB);
edist = sort(dipBS);
test.dec = dip>=edist(ihi);
test.pval = sum(edist>dip)/NB;

end

