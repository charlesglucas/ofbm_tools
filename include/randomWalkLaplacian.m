function [L,Lrw,Lsym] = randomWalkLaplacian(W)
% returns the combinatorial Laplacian L and the random walk Laplacian Lrw
% of a graph with similarity matrix W
%
% Charles-Gérard Lucas, ENS Lyon, 2021

d = sum(W,2);
D = diag(d);
L = D - W;
P = size(W,1);
Lrw = zeros(P);
Lsym = zeros(P);

for i = 1:P
    for j = 1:P
        if d(i) ~=0, Lrw(i,j) = -W(i,j)/d(i)*(i~=j) + (i==j); end
    end
end
    
for i = 1:P
    for j = 1:P
        if d(i) ~=0 && d(j) ~=0, Lsym(i,j) = -W(i,j)/(d(i)*d(j))^(1/2)*(i~=j) + (i==j); end
    end
end


end