function [nb_clusters,stats,pvals] = SilvermanModeEstim(x,alpha,params)

% returns the estimated number of clusters by iterating the Siverman's test
%
%  Input:     x:       samples 
%             alpha:   significance level of the test
%             params:  structure of parameters
%                       - tol: tolerance thresholg to compute the
%                       Silverman's statistic
%                       - nbDensity: number of points to estimate the
%                       probability density function
%
%  Output:    test:    structure of half normal testing procedures
%                       - nb_clusters: estimated number of clusters
%                       - stats: list of the computed Silverman's
%                       statistics
%                       - pvals: list of the computed Silverman's p-values
%
%
% Charles-Gérard Lucas, ENS Lyon, 2023

if ~exist('params'), params = {}; end
if ~isfield(params,'tol'), params.tol=1e-6; end
if ~isfield(params,'nbDensity'), params.nbDensity=1000; end

a = 1e-4; b = .1;
nb_clusters = 1; pvalBS = 0;
while pvalBS < alpha
    nb_clusters = nb_clusters + 1;
    stats(nb_clusters) = SilvermanStat(x,nb_clusters,[a b],params.tol);
    for r=1:100
        eps = randn(1,length(x)); zBS = (datasample(x,length(x))+stats(nb_clusters)*eps)/sqrt(1+(stats(nb_clusters)/std(x))^2);
        nbmBS(r) = sum(islocalmax(ksdensity(zBS,'NumPoints',params.nbDensity,'BandWidth',stats(nb_clusters),'function','pdf')));
    end
    pvalBS = mean(nbmBS>nb_clusters); pvals(nb_clusters) = pvalBS;
end

% compute at least 3 pvalues
if nb_clusters<3
    for nbmt = nb_clusters+1:4
        stats(nbmt) = SilvermanStat(x,nbmt,[a b],params.tol);
        for r=1:100
            eps = randn(1,length(x)); zBS = (datasample(x,length(x))+stats(nbmt)*eps)/sqrt(1+(stats(nbmt)/std(x))^2);
            nbmBS(r) = sum(islocalmax(ksdensity(zBS,'NumPoints',params.nbDensity,'BandWidth',stats(nbmt),'function','pdf')));
        end
        pvalBS = mean(nbmBS>nbmt); pvals(nbmt) = pvalBS;
    end
end

end

