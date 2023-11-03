function [dec,s,ps,index] = fdrcorrection(pvals,alpha)
% returns FDR correction
%
% Patrice Abry, Lyon, Oct 2020

m =length(pvals) ; 
[ps,index] = sort(pvals) ; 
mm = 1:m ; 
C = sum(1./mm) ;
dec = zeros(6,m);

s(1,:) = alpha*ones(size(mm)) ;
s(2,:) = alpha*mm/m ; % Benjamini-Hochberg - indep test
s(3,:) = alpha*mm/(m*C) ; % Benjamini-Yekutieli - dep test
s(4,:) = alpha/m*ones(size(mm)) ; % Bonferroni

for k=1:1:m
    dec(1,k) = ps(k)<=s(1,k) ; 
    dec(2,k) = ps(k)<=s(2,k) ; 
    dec(3,k) = ps(k)<=s(3,k) ; 
    dec(4,k) = ps(k)<=s(4,k) ; 
end

for j=2:3
tmp = dec(j,:); idx = find(tmp==1);
if ~isempty(idx), dec(j,1:idx(end)) = ones(1,idx(end)); end
end

end
