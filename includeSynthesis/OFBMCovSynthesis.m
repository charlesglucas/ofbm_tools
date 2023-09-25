% OFBM_Cov_Synthesis
% PA, June 2018, Lyon
% Design Cov for M-variate OFBM general case design

function [data,dataNoMix,Rxapprox,approxFlag,BmArray,RX,K,funch] = OFBMCovSynthesis(params) ; 
% [data,dataNoMix,params] = OFBM_Synthesis(params) ; 

name = 'Gaussian'; 
NameDemo = 'DemoGene' ; 
params.T = [0 1] ; 
data = [] ; 
dataNoMix = [] ; 

% Params Check
P0 = length(params.H) ; 
P1 = max(size(params.Correl)) ; 
P2 = max(size(params.W)) ; 
P3 = length(params.Sigma2) ; 
if ((P2~= P0)||(P1~= P0)||(P3~= P0))
    disp('/!\ Dimension Mismatch /!\');
    return
end
%
h    = params.H ;
dvect = h - 1/2 ; 
nbsamples = params.nbsamples ;  
% d = (repmat(h-0.5',1,P0) + repmat(h-0.5,P0,1))/2 ; 
Sigmae = diag(sqrt(params.Sigma2))*params.Correl*diag(sqrt(params.Sigma2)) ; 
% params.j2 = min(params.j2,log2(nbsamples) - Nwt - 3) ;
% alpha = params.alpha ;
% gamint =  0.5 ; % Ensures L2 Norm for Wav Coef 

% Compute Cov
Ncov = nbsamples+1; 

% AAt
for k1 =1:1:P0
    for k2 =1:1:P0
        M(k1,k2) = sqrt(params.Sigma2(k1)*params.Sigma2(k2))*params.Correl(k1,k2)*gamma(dvect(k1)+dvect(k2)+2)*sin(pi*(dvect(k1)+dvect(k2)+1)/2) ; 
    end
end
while det(M) <=0 
   disp('/!\ OFBM is not correctly defined /!\');
   det(M)
   break
end

%RY = CovarMfGn_sigma(Ncov,[h(1,1) h(2,2)],AAt);
RY = CovarMfGn_sigma_v2P(Ncov,h,Sigmae);

% Map RY to RX
[RX,K,funch] = MapCovar(RY,name,params.T);
% NOTE: RX is not necessarily a valid covariance structure 
[Rxapprox,approxFlag,BmArray] = InitStepMultivarGaussBestApprox(RX,nbsamples);
% [Rxapprox,approxFlag,BmArray] = InitMGBestApproxLapack(RX,nbsamples);
if approxFlag,
    disp('Approximation was needed to ensure psd of inverted target covariance.');
else
    disp('Inverted target covariance was psd - no approximation needed.');
end

P = length(RX); % number of components in the series
X = cell(1,P);

% Synthetize OfBm
for n = 1:1:params.nbcopies
for kk=1:floor(K/2),
        [Xcopy1,Xcopy2] = SynthStepMultivarGaussBestApprox(BmArray,nbsamples);        
        % Xcopy1 and Xcopy2 are independent copies of the multivariate Gaussian process 
        for p=1:P,
            X{p}(2*kk-1,:) = Xcopy1(p,:);
            X{p}(2*kk,:) = Xcopy2(p,:);        
        end
    end
    if rem(K,2)>0,
        % add one more component since K is odd
        [Xcopy1,Xcopy2] = SynthStepMultivarGaussBestApprox(BmArray,nbsamples);        
        for p=1:P,
            X{p}(K,:) = Xcopy1(p,:);
        end
    end
    % Transform Gaussian series
    for p=1:P,
        Y(p,:) = funch(X{p});
    end
     % - (5) - Mixture
      data(:,:,n) = params.W * Y ;
      dataNoMix(:,:,n) = Y ; 
      for k=1:1:P0
          data(k,:,n)  = data(k,:,n) * sqrt(params.Sigma2(k)) ;
      end
      
end
