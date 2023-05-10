function [W_corr,decL] = SingleExpDetection(W,alpha,corr)
% returns weight matrix after detection of single exponent
%
%  Input:     est:     structure of estimation method 
%             alpha:   significance level of the tests 
%             corr:    correction procedure
%                           1 -  No correction
%                           2 -  Benjamini-Hochberg
%                           3 -  Benjamini-Yekutieli
%                           4 -  Bonferroni
%
%  Output:    W_corr:  similarity matrix
%             decL:    pairwise corrected decisions
%
% Charles-Gérard Lucas, ENS Lyon, 2021

P = size(W,1);
k = 0 ;
for k1=1:P-1
    for k2=k1+1:P
        k = k+1;
        pvals(k) = W(k1,k2);
    end
end
    
Adj2 = zeros(P);

[dec,~,~,index] = fdrcorrection(pvals,alpha);
decs(index) = dec(corr,:); 

k = 0;
for k1=1:P-1
    for k2=k1+1:P
        k = k+1;
        decL(k1,k2) = decs(k);
        Adj2(k1,k2) = 1-decL(k1,k2);
        Adj2(k2,k1) = 1-decL(k1,k2);
    end
end

nb_iso = sum((sum(Adj2,2)==0));

% corrected adjacency
d = sum(W,2);  [~,id] = sort(d);
W_corr = W;
W_corr(:,id(1:nb_iso)) = zeros(P,nb_iso);
W_corr(id(1:nb_iso),:) = zeros(nb_iso,P);

end

