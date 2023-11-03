function [nbcluster, cluster, eigenvalues] = autoClustering(L)
% returns the number of clusters and clusters using spectral clustering
% on a graph with Laplacian matrix L
%
% Charles-Gérard Lucas, ENS Lyon, 2021

% compute eigenvalues and eigenvactors of the Laplacian
[V,DV] = eig(L); 
[eigenvalues,id] = sort(diag(DV));
vec = V(:,id);

% number of clusters with maximum eigengap
if eigenvalues(end) ~= 0
    [~,nbcluster] = max(diff(eigenvalues));
else 
    nbcluster = size(L,1);
end
% k-means
jc = nbcluster;
cluster = kmeans(real(vec(:,1:jc)),jc)';

end