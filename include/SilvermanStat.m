function [hcrit] = SilvermanStat(samples,nb_mode,range,tol)

% returns the Silverman's statistic
%
% Charles-Gérard Lucas, ENS Lyon, 2023

a = range(1); b = range(2);
while b-a > tol
    h = (b+a)/2;
    nbm = sum(islocalmax(ksdensity(samples,0:1/1000:1,'BandWidth',h,'function','pdf')));
    if nbm > nb_mode, a = h; else, b = h; end
end
hcrit = b;

end

