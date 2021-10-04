function [WD] = Dx_bootstrap_OFBM(WD,NB,LB,j1,IID)

%% If you use this code, please cite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H. Wendt, P. Abry, G. Didier, "Wavelet domain bootstrap for testing the
% equality of bivariate self-similarity exponents,"	IEEE Workshop
% Statistical Signal Proces. (SSP), Freiburg, Germany, June 2018.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3; LB=1; end
if nargin<4; j1=1; end
if nargin<5; IID=0; end

%if ~IID; rng('shuffle'); rst=rng; end
if ~IID; rst=rng; end

J=length(WD{1});
for k=1:length(WD);
    % use same draws for each component
    if ~IID; rng(rst); end
    for j=j1:J
        tmp=WD{k}(j).value_noabs;
        nj=length(tmp);
        %ii=ceil(rand(NB,nj)*nj);
        if LB==1
        	ii0=ceil(rand(1,NB*nj)*nj); ii=reshape(ii0,NB,nj);
        else
            LBj=LB; if LBj>nj/2; LBj=round(nj/4); end
            NBB=max(4,floor(nj/LBj));
            ii0=ceil(rand(1,NB*NBB)*(nj-LBj));
            ii1=repmat(ii0,LBj,1)+repmat((1:LBj)',1,NB*NBB);
            ii=reshape(ii1,NB,[]);
        end
        try
        WD{k}(j).value_noabs_bs=WD{k}(j).value_noabs(ii);
        catch
            keyboard
        end
    end
end
