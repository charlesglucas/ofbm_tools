function [test] = BSFoldedNormalTest(est,alpha)

% returns the half normal pairwise testing procedure with parameter
% estimating with a bootstrap procedure reproducing the observed statistic
%
%  Input:     est:     structure of estimation method 
%             alpha:   significance level of the tests  
%
%  Output:    test:    structure of folded normal testing procedures
%                       - SBE: estimated test parameter
%                       - pval: pairwise p-values
%                       - decpw: pairwise decisions
%                       - decBFpw: pairwise decisions corrected with
%                       Bonferroni procedure
%                       - decHocpw: pariwise decisions corrected with
%                       Benjamini-Hochberg procedure
%                       - decYekpw: pairwise decisions corrected with
%                       Benjamini-Yekutieli procedure
%                       - decNC: decision without correction
%                       - decBF: decision from Bonferroni correction
%                       - decHoc: decision from Benjamini-Hochberg
%                       correction
%                       - decYek: decision from Benjamini-Yekutieli
%                       correction
%
%
% Charles-Gérard Lucas, ENS Lyon, 2021

P = size(est.lambdaBS,2);
tc = icdf('HalfNormal',1-alpha,0,1);

lambdaBSsort=sort(est.lambdaBS,2,'ascend');

% pairwise decisions and p-value
for p1=1:1:P-1
    for p2=p1+1:P
        [~,sigma] = FNparameter(lambdaBSsort(:,p2)-lambdaBSsort(:,p1));
        test.SBE(p1,p2) = sigma;
        pvals(p1,p2)=1-cdf('HalfNormal',abs(est.lambdasort(p1)-est.lambdasort(p2))/test.SBE(p1,p2),0,1);
    end
end

% FDR correction
PvalSeq = pvals;
PvalSeq(P,:) = zeros(1,P);
test.pval = diag(PvalSeq,1);
[dec,~,~,index] = fdrcorrection(test.pval,alpha);
for r=1:6, for k=1:P-1, decs(r,index(k)) = dec(r,k); end; end
test.decpw = decs(1,:);
test.decBFpw = decs(4,:);
test.decHocpw = decs(2,:);
test.decYekpw = decs(3,:);
test.decNC = pairwiseTestProduct(decs(1,:));
test.decBF = pairwiseTestProduct(decs(4,:));
test.decHoc = pairwiseTestProduct(decs(2,:));
test.decYek = pairwiseTestProduct(decs(3,:));

end

