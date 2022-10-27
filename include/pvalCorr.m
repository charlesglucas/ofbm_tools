function [Adj_corr,decL] = pvalCorr(Adj,alpha,r)
% returns weight matrix after detection of single exponent
%
% Charles-Gérard Lucas, ENS Lyon, 2021

P = size(Adj,1);
k = 0 ;
for k1=1:P-1
    for k2=k1+1:P
        k = k+1;
        pvals(k) = Adj(k1,k2);
    end
end
    
Adj2 = zeros(P);

[dec,~,~,index] = fdrcorrection(pvals,alpha);
decs(index) = dec(r,:); 

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
d = sum(Adj,2);  [~,id] = sort(d);
Adj_corr = Adj;
Adj_corr(:,id(1:nb_iso)) = zeros(P,nb_iso);
Adj_corr(id(1:nb_iso),:) = zeros(nb_iso,P);

end

