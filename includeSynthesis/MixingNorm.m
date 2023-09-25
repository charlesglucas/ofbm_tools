
function W = MixingNorm(P) ; 

W=(P) ; 
for i=1:length(P); 
    W(i,i)= 1 ; 
    W(:,i)=W(:,i)/sqrt(sum((W(:,i)).^2)) ; 
end; 


