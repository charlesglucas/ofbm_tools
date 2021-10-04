function [decSingle] = pairwiseTestProduct(dec)
% rejects the null hypothesis H_1 = ... = H_M if any of the pairwise 
% decisions rejects the null hypothesis H_i = H_j

dec0 = 1;
for i=1:length(dec)
    dec0 = dec0*(1-dec(i));
end
decSingle = 1-dec0;

end
