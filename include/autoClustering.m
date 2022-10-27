function [nbcluster, cluster, eigenvalues] = autoClustering(L)
% returns the number of clusters and clusters using spectral clustering
% on a graph with Laplacian matrix L
%
% Charles-Gérard Lucas, ENS Lyon, 2021

[V,DV] = eig(L); 
% number of clusters with maximum eigengap
[eigenvalues,id] = sort(diag(DV));
vec = V(:,id);
if eigenvalues(end) ~= 0
    [~,nbcluster] = max(diff(eigenvalues));
else 
    nbcluster = size(L,1);
end
% k-means
jc = nbcluster;
cluster = kmeans(real(vec(:,1:jc)),jc)';

end