function [estT] = OFBM_estimBC_BS_test(est,alpha)
%  Input:     est:     structure of estimation method 
%             alpha:   significance level of the tests   
%
%  Output:    estT:    structure of estimation method with testing
%                   procedures in addition
%
% Charles-Gérard Lucas, ENS Lyon, 2021

P = size(est.lambdaBS,2);
NB = size(est.lambdaBS,1);
tc=norminv(1-alpha/2);
tc2 = icdf('HalfNormal',1-alpha,0,1);
ilo=round(alpha/2*NB);
ihi=round((1-alpha/2)*NB);
ihi2=round((1-alpha)*NB);

estT.lambdaBSsort=sort(est.lambdaBS-mean(est.lambdaBS),2,'ascend');
estT.lambdaBSsortH1=sort(est.lambdaBS,2,'ascend');

p = 0;
for p1=1:1:P-1
    for p2=p1+1:P
        % bootstrap variances of delta
        estT.VarDelta(p1,p2)=var(est.lambdaBS(:,p1)-est.lambdaBS(:,p2));
        estT.VarDeltasort(p1,p2)=var(estT.lambdaBSsort(:,p1)-estT.lambdaBSsort(:,p2));
        estT.MeanDelta(p1,p2)=mean(est.lambdaBS(:,p1)-est.lambdaBS(:,p2));
        estT.SBE(p1,p2)=sqrt(estT.VarDelta(p1,p2));
        estT.SBEsort(p1,p2) = sqrt(estT.VarDeltasort(p1,p2))/sqrt(1-2/pi);
        [~,sigma] = FNparameter(estT.lambdaBSsortH1(:,p2)-estT.lambdaBSsortH1(:,p1));
        estT.SBEsort_v2(p1,p2) = sigma;
        
        % Gaussian and half-Gaussian tests
        estT.dec(p1,p2)=abs(est.lambda(p1)-est.lambda(p2))>tc*estT.SBE(p1,p2);
        estT.pval(p1,p2)=2*(1-normcdf(abs(est.lambda(p1)-est.lambda(p2))/estT.SBE(p1,p2)));
        estT.decsort(p1,p2)=abs(est.lambdasort(p1)-est.lambdasort(p2))>tc2*estT.SBEsort(p1,p2);
        estT.pvalsort(p1,p2)=1-cdf('HalfNormal',abs(est.lambdasort(p1)-est.lambdasort(p2))/estT.SBEsort(p1,p2),0,1);
        estT.decsort_v2(p1,p2)=abs(est.lambdasort(p1)-est.lambdasort(p2))>tc2*estT.SBEsort_v2(p1,p2);
        estT.pvalsort_v2(p1,p2)=1-cdf('HalfNormal',abs(est.lambdasort(p1)-est.lambdasort(p2))/estT.SBEsort_v2(p1,p2),0,1);
        
        % percentile Bootstrap test 
        eD=(est.lambda(p1)-est.lambda(p2));
        edist=sort(est.lambdaBS(:,p1)-est.lambdaBS(:,p2) - mean(est.lambdaBS(:,p1)-est.lambdaBS(:,p2)));
        estT.decBS(p1,p2)= eD>=edist(ihi) | eD<=edist(ilo);
        QT=sum(edist<eD);
        estT.pvalBS(p1,p2)=2/NB*min(QT,NB-QT);
        
        % ordered percentile Bootstrap test 
        eDsort=(est.lambdasort(p2)-est.lambdasort(p1));
        edistsort=sort(estT.lambdaBSsort(:,p2)-estT.lambdaBSsort(:,p1));
        estT.decBSsort(p1,p2)= eDsort>=edistsort(ihi2);
        QT=sum(edistsort<eDsort);
        estT.pvalBSsort(p1,p2)=1/NB*(NB-QT);
        
        invSigma = inv(cov(est.lambdaBS(:,p1),est.lambdaBS(:,p2)));
        cL = est.lambda([p1,p2]) - mean(est.lambda([p1,p2]));
        estT.Tpw(p1,p2) = sum(cL*invSigma*cL');
        estT.decTpw(p1,p2) = estT.Tpw(p1,p2) >= chi2inv(1-alpha,1);
        estT.pvalTpw(p1,p2) = 1 - chi2cdf(estT.Tpw(p1,p2),1);
    end
end

% FDR correction for half-normal test
PvalSeq = estT.pvalsort;
PvalSeq(P,:) = zeros(1,P);
[dec,~,~,index] = fdrcorrection(diag(PvalSeq,1),alpha);
for r=1:6, for k=1:P-1, decs(r,index(k)) = dec(r,k); end; end
estT.decsortpw = decs(1,:);
estT.decsortBFpw = decs(4,:);
estT.decsortHocpw = decs(2,:);
estT.decsortYekpw = decs(3,:);
estT.decsortHoc2pw = decs(5,:);
estT.decsortYek2pw = decs(6,:);
estT.decsort2 = pairwiseTestProduct(decs(1,:));
estT.decsortBF = pairwiseTestProduct(decs(4,:));
estT.decsortHoc = pairwiseTestProduct(decs(2,:));
estT.decsortYek = pairwiseTestProduct(decs(3,:));
estT.decsortHoc2 = pairwiseTestProduct(decs(5,:));
estT.decsortYek2 = pairwiseTestProduct(decs(6,:));

% FDR correction for alternative half-normal test
% using folded normal parameter estimation
PvalSeq = estT.pvalsort_v2;
PvalSeq(P,:) = zeros(1,P);
[dec,~,~,index] = fdrcorrection(diag(PvalSeq,1),alpha);
for r=1:6, for k=1:P-1, decs(r,index(k)) = dec(r,k); end; end
estT.decsortpw_v2 = decs(1,:);
estT.decsortBFpw_v2 = decs(4,:);
estT.decsortHocpw_v2 = decs(2,:);
estT.decsortYekpw_v2 = decs(3,:);
estT.decsortHoc2pw_v2 = decs(5,:);
estT.decsortYek2pw_v2 = decs(6,:);
estT.decsort2_v2 = pairwiseTestProduct(decs(1,:));
estT.decsortBF_v2 = pairwiseTestProduct(decs(4,:));
estT.decsortHoc_v2 = pairwiseTestProduct(decs(2,:));
estT.decsortYek_v2 = pairwiseTestProduct(decs(3,:));
estT.decsortHoc2_v2 = pairwiseTestProduct(decs(5,:));
estT.decsortYek2_v2 = pairwiseTestProduct(decs(6,:));

% FDR correction for non parametric sorted test
PvalSeq = estT.pvalBSsort;
PvalSeq(P,:) = zeros(1,P);
[dec,~,~,index] = fdrcorrection( diag(PvalSeq,1),alpha);
for r=1:4, for k=1:P-1, decs(r,index(k)) = dec(r,k); end; end
estT.decsortBSpw = decs(1,:);
estT.decsortBFBSpw = decs(4,:);
estT.decsortHocBSpw = decs(2,:);
estT.decsortYekBSpw = decs(3,:);
estT.decsortHoc2BSpw = decs(5,:);
estT.decsortYek2BSpw = decs(6,:);
estT.decsortBS2 = pairwiseTestProduct(decs(1,:));
estT.decsortBFBS = pairwiseTestProduct(decs(4,:));
estT.decsortHocBS = pairwiseTestProduct(decs(2,:));
estT.decsortYekBS = pairwiseTestProduct(decs(3,:));
estT.decsortHoc2BS = pairwiseTestProduct(decs(5,:));
estT.decsortYek2BS = pairwiseTestProduct(decs(6,:));

% FDR correction for normal test
PvalSeq = estT.pval;
PvalSeq(P,:) = zeros(1,P);
[dec,~,~,index] = fdrcorrection(diag(PvalSeq,1),alpha);
for r=1:4, for k=1:P-1, decs(r,index(k)) = dec(r,k); end; end
estT.decpw = decs(1,:);
estT.decBFpw = decs(4,:);
estT.decHocpw = decs(2,:);
estT.decYekpw = decs(3,:);
estT.decHoc2pw = decs(5,:);
estT.dec2 = pairwiseTestProduct(decs(1,:));
estT.decBF = pairwiseTestProduct(decs(4,:));
estT.decHoc = pairwiseTestProduct(decs(2,:));
estT.decYek = pairwiseTestProduct(decs(3,:));
estT.decHoc2 = pairwiseTestProduct(decs(5,:));
estT.decYek2 = pairwiseTestProduct(decs(6,:));

% FDR correction for non parametric test
PvalSeq = estT.pvalBS;
PvalSeq(P,:) = zeros(1,P);
estT.pvalseq = diag(PvalSeq,1);
[dec,~,~,index] = fdrcorrection(estT.pvalseq,alpha);
for r=1:4, for k=1:P-1, decs(r,index(k)) = dec(r,k); end; end
estT.decBSpw = decs(1,:);
estT.decBFBSpw = decs(4,:);
estT.decHocBSpw = decs(2,:);
estT.decYekBSpw = decs(3,:);
estT.decHoc2BSpw = decs(5,:);
estT.decBS2 = pairwiseTestProduct(decs(1,:));
estT.decBFBS = pairwiseTestProduct(decs(4,:));
estT.decHocBS = pairwiseTestProduct(decs(2,:));
estT.decYekBS = pairwiseTestProduct(decs(3,:));
estT.decHoc2BS = pairwiseTestProduct(decs(5,:));
estT.decYek2BS = pairwiseTestProduct(decs(6,:));

% compute approximate chi2 test
estT.invSigma = inv(cov(est.lambdaBS));
cL = est.lambda - mean(est.lambda);
estT.T=cL*estT.invSigma*cL';
estT.T2=sum( (est.lambda-mean(est.lambda)).^2 );
estT.Tdec = estT.T >= chi2inv(1-alpha,P-1);
estT.Tpval = 1 - chi2cdf(estT.T,P-1);

end

