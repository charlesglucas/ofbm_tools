function [estT] = OFBM_estimBC_BS_test(est,aa,params)
%  Input:     est:     structure of estimation method 
%             alpha:   significance level of the tests     
%             params:  structure of parameters
%
%  Output:    estT:    structure of estimation method with testing
%                   procedures in adddition
%
% Charles-Gérard Lucas, ENS Lyon, 2021

estT = est;
P = params.P;
NB = params.NB;
tc=norminv(1-aa/2);
%tcbis = icdf('normal',1-aa/2,0,1);
%tc2=norminv(1-aa);
tc2 = icdf('HalfNormal',1-aa,0,1);
ilo=round(aa/2*NB);
ihi=round((1-aa/2)*NB);
ihi2=round((1-aa)*NB);

estT.lambdaBSsort=sort(estT.lambdaBS-mean(estT.lambdaBS),2,'ascend');
%estT.lambdaBSsort=sort(estT.lambdaBS,2,'ascend');

p = 0;
for p1=1:1:P-1
    for p2=p1+1:P
        % bootstrap variances of delta
        estT.VarDelta(p1,p2)=var(estT.lambdaBS(:,p1)-estT.lambdaBS(:,p2));
        estT.VarDeltasort(p1,p2)=var(estT.lambdaBSsort(:,p1)-estT.lambdaBSsort(:,p2));
        estT.MeanDelta(p1,p2)=mean(estT.lambdaBS(:,p1)-estT.lambdaBS(:,p2));
        estT.SBE(p1,p2)=sqrt(estT.VarDelta(p1,p2));
        estT.SBEsort(p1,p2) = sqrt(estT.VarDeltasort(p1,p2))/sqrt(1-2/pi);
        
        % Gaussian and half-Gaussian tests
        estT.dec(p1,p2)=abs(estT.lambda(p1)-estT.lambda(p2))>tc*estT.SBE(p1,p2);
        estT.pval(p1,p2)=2*(1-normcdf(abs(estT.lambda(p1)-estT.lambda(p2))/estT.SBE(p1,p2)));
        %estT.pvalbis(p1,p2)=2*(1-cdf('Normal',abs(estT.lambda(p1)-estT.lambda(p2))/estT.SBE(p1,p2),0,1));
        estT.decsort(p1,p2)=abs(estT.lambdasort(p1)-estT.lambdasort(p2))>tc2*estT.SBEsort(p1,p2);
        %estT.decsortbis(p1,p2)=abs(estT.lambdasort(p1)-estT.lambdasort(p2))>tc2*estT.SBE(p1,p2);
        estT.pvalsort(p1,p2)=1-cdf('HalfNormal',abs(estT.lambdasort(p1)-estT.lambdasort(p2))/estT.SBEsort(p1,p2),0,1);
        %estT.pvalsortbis(p1,p2)=1-cdf('HalfNormal'ù,abs(estT.lambdasort(p1)-estT.lambdasort(p2))/estT.SBE(p1,p2),0,1);
        
        % percentile Bootstrap test 
        eD=(estT.lambda(p1)-estT.lambda(p2));
        edist=sort(estT.lambdaBS(:,p1)-estT.lambdaBS(:,p2) - mean(estT.lambdaBS(:,p1)-estT.lambdaBS(:,p2)));
        estT.decBS(p1,p2)= eD>=edist(ihi) | eD<=edist(ilo);
        QT=sum(edist<eD);
        estT.pvalBS(p1,p2)=2/NB*min(QT,NB-QT);
        
        % ordered percentile Bootstrap test 
        eDsort=(estT.lambdasort(p2)-estT.lambdasort(p1));
        edistsort=sort(estT.lambdaBSsort(:,p2)-estT.lambdaBSsort(:,p1));
        estT.decBSsort(p1,p2)= eDsort>=edistsort(ihi2);
        QT=sum(edistsort<eDsort);
        estT.pvalBSsort(p1,p2)=1/NB*(NB-QT);
        
        % ordered percentile Bootstrap test with permutation
%         edistsortbis=sort((estT.lambdaBS(:,estT.sortid(p2))-estT.lambdaBS(:,estT.sortid(p1))) - mean(estT.lambdaBS(:,estT.sortid(p2))-estT.lambdaBS(:,estT.sortid(p1))));
%         estT.decBSsortbis(p1,p2)= eDsort>=edistsortbis(ihi2);
%         QT=sum(edistsortbis<eDsort);
%         estT.pvalBSsortbis(p1,p2)=1/NB*(NB-QT);
        
        invSigma = inv(cov(estT.lambdaBS(:,p1),estT.lambdaBS(:,p2)));
        cL = estT.lambda([p1,p2]) - mean(estT.lambda([p1,p2]));
        estT.Tpw(p1,p2) = sum(cL*invSigma*cL');
        estT.decTpw(p1,p2) = estT.Tpw(p1,p2) >= chi2inv(1-aa,1);
        estT.pvalTpw(p1,p2) = 1 - chi2cdf(estT.Tpw(p1,p2),1);
    end
end

% FDR correction for half-normal test
PvalSeq = estT.pvalsort;
PvalSeq(P,:) = zeros(1,P);
[dec,~,~,index] = fdrcorrection(diag(PvalSeq,1),aa);
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

% FDR correction for non parametric sorted test
PvalSeq = estT.pvalBSsort;
PvalSeq(P,:) = zeros(1,P);
[dec,~,~,index] = fdrcorrection( diag(PvalSeq,1),aa);
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
[dec,~,~,index] = fdrcorrection(diag(PvalSeq,1),aa);
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
[dec,~,~,index] = fdrcorrection(estT.pvalseq,aa);
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

% extract test h1<h2<...<hM 
% tmp=estT.decsort; estT.decdiag=diag(tmp(:,2:end));
% tmp=estT.pvalsort; estT.pvaldiag=diag(tmp(:,2:end));
% tmp=estT.decBSsort; estT.decBSdiag=diag(tmp(:,2:end));
% tmp=estT.pvalBSsort; estT.pvalBSdiag=diag(tmp(:,2:end));
% tmp=estT.decBSsortbis; estT.decBSdiagbis=diag(tmp(:,2:end));
% tmp=estT.pvalBSsortbis; estT.pvalBSdiagbis=diag(tmp(:,2:end));

% compute approximate chi2 test
estT.invSigma = inv(cov(estT.lambdaBS));
cL = estT.lambda - mean(estT.lambda);
estT.T=cL*estT.invSigma*cL';
estT.T2=sum( (estT.lambda-mean(estT.lambda)).^2 );
estT.Tdec = estT.T >= chi2inv(1-aa,P-1);
estT.Tpval = 1 - chi2cdf(estT.T,P-1);

% Chi2 test with sigma diagonal
% estT.Tbis=sum( ( estT.lambda - mean(estT.lambda) ).^2 ./ estT.varBS );
% estT.Tbisdec = estT.Tbis >= chi2inv(1-aa,P-1);
% estT.Tbispval = 1 - chi2cdf(estT.Tbis,P-1);

% for p1=1:1:P
%     studBS(:,p1)=(est.lambdaBS(:,p1)-mean(est.lambdaBS(:,p1)))./std(est.lambdaBS(:,p1));
%     studBSbis(:,p1)=studBS(:,p1)+mean(est.lambdaBS(:,p1))-mean(est.lambda);
% end
% aa=norminv(0.75);

% t test ?
% for p1=1:1:P
%     %elBS(:,p1)=est.lambdaBS(:,p1)-mean(est.lambdaBS(:,p1));
%     elBS(:,p1)=estT.lambdaBS(:,p1)-estT.lambda(p1);
% end
% mBS=mean(elBS,2);
% tBS=zeros(size(mBS));
% for p1=1:1:P
%     %tBS=tBS+elBS(:,p1).^2;
%     tBS=tBS+(elBS(:,p1)-mBS).^2;
% end
% estT.tBS=sort(tBS);
% estT.t=sum((estT.lambda-mean(estT.lambda)).^2);
% estT.tdec= estT.t>=estT.tBS(ihi2);
% estT.tpval=sum(estT.tBS>estT.t)/NB;

end

