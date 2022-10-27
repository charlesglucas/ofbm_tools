function [avH] = averagedClusters(H,cluster)
% returns self-similarity exponents averaged in each cluster
%
% Charles-Gérard Lucas, ENS Lyon, 2021

avH = H;
for i=1:max(cluster)
    idx = find(cluster == i);
    avH(idx) = mean(avH(idx));
end

end

