% OFBMSynthesis.m
% PA, Lyon, July 2019

function [data,dataNoMix] = OFBMSynthesis(K,RX,BmArray,funch,params) ; 

%[~,~,Rxapprox,approxFlag,BmArray,RX,K,funch] = OFBMCovSynthesis(params) ; 

P = params.P ; 
if params.CorrDiag == 0
    P = length(RX); % number of components in the series
    X = cell(1,P);
    for kk=1:floor(K/2),
        [Xcopy1,Xcopy2] = SynthStepMultivarGaussBestApprox(BmArray,params.nbsamples);
        % Xcopy1 and Xcopy2 are independent copies of the multivariate Gaussian process
        for p=1:P,
            X{p}(2*kk-1,:) = Xcopy1(p,:);
            X{p}(2*kk,:) = Xcopy2(p,:);
        end
    end
    if rem(K,2)>0,
        % add one more component since K is odd
        [Xcopy1,Xcopy2] = SynthStepMultivarGaussBestApprox(BmArray,params.nbsamples);
        for p=1:P,
            X{p}(K,:) = Xcopy1(p,:);
        end
    end
    % Transform Gaussian series
    for p=1:P,
        Y(p,:) = funch(X{p});
    end
else
    for p=1:P,
        [tmp1,Y(p,:),tmp2] = synthfbmcircul(params.nbsamples,params.H(p)) ; 
    end
    
end
for k=1:1:P
    Y(k,:)  = Y(k,:)./std(Y(k,:))*sqrt(params.Sigma2(k)) ;
end
% - (5) - Mixture
data = params.W * Y ;
dataNoMix = Y ; 
if params.FBM == 1
    for k=1:1:P
        data(k,:)  = cumsum(data(k,:)) ;
        dataNoMix(k,:) = cumsum(dataNoMix(k,:)) ;
    end
end




