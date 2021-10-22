function [WD,WDD] = Dx_bootstrap_OFBM(WD,NB,LB,j1,IID)

dbs=0; if nargout>1; dbs=1; end

if nargin<3; LB=1; end
if nargin<4; j1=1; end
if nargin<5; IID=0; end

%if ~IID; rng('shuffle'); rst=rng; end
if ~IID; rst=rng; end

WDD=[];
J=length(WD{1});
for k=1:length(WD);
    % use same draws for each component
    if ~IID; rng(rst); end
    for j=j1:J
%         try
            tmp=WD{k}(j).value_noabs;
%         catch
%             tmp=WD{k}(j).value;
%         end
        nj=length(tmp);
        %ii=ceil(rand(NB,nj)*nj);
        if LB==1
        	ii0=ceil(rand(1,NB*nj)*nj); ii=reshape(ii0,NB,nj);
        else
            LBj=LB; if LBj>nj/2; LBj=round(nj/4); end
            %NBB=max(4,floor(nj/LBj));
            NBB=max(4,ceil(nj/LBj));
            ii0=ceil(rand(1,NB*NBB)*(nj-LBj));
            ii1=repmat(ii0,LBj,1)+repmat((1:LBj)',1,NB*NBB);
            %ii=reshape(ii1,NB,[]);
            ii=reshape(ii1,[],NB)';
            if dbs>0 % do single sample dbs-fold bootstrap
                madds=repmat((0:NB-1)*NBB,NBB,1);
                vadds=madds(:)';
                ii0s=ii0;
                for idb=1:dbs
                    subi=ceil(rand(1,NB*NBB)*NBB)+vadds;
                    ii0s=ii0s(subi);
                end
                ii1s=repmat(ii0s,LBj,1)+repmat((1:LBj)',1,NB*NBB);
                iis=reshape(ii1s,[],NB)';
            end
        end
        try
            WD{k}(j).value_noabs_bs=WD{k}(j).value_noabs(ii);
            if dbs==1
                WD{k}(j).value_noabs_doublebs=WD{k}(j).value_noabs(iis);
            end
        catch
            keyboard
        end
    end
end
