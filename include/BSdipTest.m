function [dec] = BSdipTest(estbc,alpha)

NB = size(estbc.lambdaBS,1);
dip = HartigansDipTest(sort(estbc.h));
for r=1:NB
	dipBS(r) = HartigansDipTest(sort(estbc.lambdaBS(r,:)-mean(estbc.lambdaBS(r,:))));
end
ihi=round((1-alpha)*NB);
edist=sort(dipBS);
dec = dip>=edist(ihi);

end

