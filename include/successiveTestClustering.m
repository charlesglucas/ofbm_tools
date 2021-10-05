function [nbcluster, cluster] = successiveTestClustering(dec)
% Charles-Gérard Lucas, ENS Lyon, August 2021

K = length(dec);
decList = dec;
cluster = zeros(1,K); cluster(1) = 1;
i = 1;
for k = 1:K
    if decList(k) ==  0
        cluster(k+1) = i;
    else
        i = i+1;
        cluster(k+1) = i;
    end
end
nbcluster = i;

end