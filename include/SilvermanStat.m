function [hcrit,nbm] = SilvermanStat(samples,nb_mode,range,tol)

a = range(1); b = range(2);
while b-a > tol
    h = (b+a)/2;
    nbm = sum(islocalmax(ksdensity(samples,0:1/1000:1,'BandWidth',h,'function','pdf')));
    if nbm > nb_mode, a = h; else, b = h; end
end
hcrit = b;

end

