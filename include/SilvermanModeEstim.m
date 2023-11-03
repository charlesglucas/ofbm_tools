function [nbm,s,pval] = SilvermanModeEstim(x,alpha,params)

if ~exist('params'), params = {}; end
if ~isfield(params,'tol'), params.tol=1e-6; end
if ~isfield(params,'nbDensity'), params.nbDensity=1000; end

a = 1e-4; b = .1;
nbm = 1; pvalBS = 0;
while pvalBS < alpha
    nbm = nbm + 1;
    s(nbm) = SilvermanStat(x,nbm,[a b],params.tol);
    for r=1:100
        eps = randn(1,length(x)); zBS = (datasample(x,length(x))+s(nbm)*eps)/sqrt(1+(s(nbm)/std(x))^2);
        nbmBS(r) = sum(islocalmax(ksdensity(zBS,'NumPoints',params.nbDensity,'BandWidth',s(nbm),'function','pdf')));
    end
    pvalBS = mean(nbmBS>nbm); pval(nbm) = pvalBS;
end

% compute at least 3 pvalues
if nbm<3
    for nbmt = nbm+1:4
        s(nbmt) = SilvermanStat(x,nbmt,[a b],params.tol);
        for r=1:100
            eps = randn(1,length(x)); zBS = (datasample(x,length(x))+s(nbmt)*eps)/sqrt(1+(s(nbmt)/std(x))^2);
            nbmBS(r) = sum(islocalmax(ksdensity(zBS,'NumPoints',params.nbDensity,'BandWidth',s(nbmt),'function','pdf')));
        end
        pvalBS = mean(nbmBS>nbmt); pval(nbmt) = pvalBS;
    end
end

end

