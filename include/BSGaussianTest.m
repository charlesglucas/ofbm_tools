function [test] = BSGaussianTest(est,alpha)
%  Input:     est:     structure of estimation method 
%             alpha:   significance level of the tests  
%
%  Output:    test:    structure of half normal testing procedures
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

P = length(est.h);
NB = size(est.lambdaBS,1);
tc=norminv(1-alpha/2);

for p1=1:1:P-1
    for p2=p1+1:P
        test.VarDelta(p1,p2)=var(est.lambdaBS(:,p1)-est.lambdaBS(:,p2));
        test.SBE(p1,p2)=sqrt(test.VarDelta(p1,p2));
        test.dec(p1,p2)=abs(est.lambda(p1)-est.lambda(p2))>tc*test.SBE(p1,p2);
        test.pval(p1,p2)=2*(1-normcdf(abs(est.lambda(p1)-est.lambda(p2))/test.SBE(p1,p2)));
    end
end

% FDR correction for normal test
k = 0; for k1=1:P-1, for k2=k1+1:P, k = k+1; PvalSeq2(k) = test.pval(k1,k2); end; end
[dec,~,~,index] = fdrcorrection(PvalSeq2,alpha);
for r=1:6, for k=1:P*(P-1)/2, decs2(r,index(k)) = dec(r,k); end; end
k = 0;
for k1=1:P-1
    for k2=k1+1:P
        k = k+1;
        test.decpw(k1,k2) = decs2(1,k);
        test.decBFpw(k1,k2) = decs2(4,k);
        test.decHocpw(k1,k2) = decs2(2,k);
        test.decYekpw(k1,k2) = decs2(3,k);
    end
end
test.dec2 = pairwiseTestProduct(decs2(1,:));
test.decBF = pairwiseTestProduct(decs2(4,:));
test.decHoc = pairwiseTestProduct(decs2(2,:));
test.decYek = pairwiseTestProduct(decs2(3,:));

