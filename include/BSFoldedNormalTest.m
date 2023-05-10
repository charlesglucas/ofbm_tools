function [test] = BSFoldedNormalTest(est,alpha)

% returns the half normal pairwise testing procedure with parameter
% estimating with a bootstrap procedure reproducing the observed statistic
%
%  Input:     est:     structure of estimation method 
%             alpha:   significance level of the tests  
%
%  Output:    test:    structure of folded normal testing procedures
%                       - decsort : pairwise decisions
%                       - pvalsort : pariwise p-values
%                       - decsortBFpw : pairwise decisions corrected with
%                       Bonferroni procedures
%                       - decsortHocpw : pariwise p-values corrected with
%                       Benjamini-Hochberg procedures
%                       - decsortYekpw : pairwise decisions corrected with
%                       Benjamini-Yekutieli procedures
%
% Charles-Gérard Lucas, ENS Lyon, 2021

P = size(est.lambdaBS,2);
tc = icdf('HalfNormal',1-alpha,0,1);

test.lambdaBSsort=sort(est.lambdaBS,2,'ascend');

% pairwise decisions and p-value
for p1=1:1:P-1
    for p2=p1+1:P
        [~,sigma] = FNparameter(test.lambdaBSsort(:,p2)-test.lambdaBSsort(:,p1));
        test.SBEsort(p1,p2) = sigma;
        test.decsort(p1,p2)=abs(est.lambdasort(p1)-est.lambdasort(p2))>tc*test.SBEsort(p1,p2);
        test.pvalsort(p1,p2)=1-cdf('HalfNormal',abs(est.lambdasort(p1)-est.lambdasort(p2))/test.SBEsort(p1,p2),0,1);
    end
end

% FDR correction
PvalSeq = test.pvalsort;
PvalSeq(P,:) = zeros(1,P);
[dec,~,~,index] = fdrcorrection(diag(PvalSeq,1),alpha);
for r=1:6, for k=1:P-1, decs(r,index(k)) = dec(r,k); end; end
test.decsortpw = decs(1,:);
test.decsortBFpw = decs(4,:);
test.decsortHocpw = decs(2,:);
test.decsortYekpw = decs(3,:);
test.decsortNC = pairwiseTestProduct(decs(1,:));
test.decsortBF = pairwiseTestProduct(decs(4,:));
test.decsortHoc = pairwiseTestProduct(decs(2,:));
test.decsortYek = pairwiseTestProduct(decs(3,:));

end

