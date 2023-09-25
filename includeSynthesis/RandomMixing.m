
function W = RandomMixing(P,sigma,seed) ; 

% sigma in (0,1) 
if nargin < 2
    sigma = 1 ; 
end
try 
rng(seed);
catch 
end
W = (rand(P,P)-0.5)*2*sigma ; 
for i=1:P
    W(i,i)= 1 ; 
%    W(:,i)=W(:,i)/sum(abs(W(:,i))); 
end; 


