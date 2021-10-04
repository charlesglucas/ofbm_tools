function [slope,intercept] = SmartLDplot(yj,nj,varyj,j1,j2,wtype)
%
%  Input:     
%             yj:       vector or matrix of structure functions
%                       if matrix: octaves j must be last dimension of matrix!!
%             varyj:    The variances of yj, precalculated (for wtype=2, otherwise ignored). 
%             nj:       the vector[1..scalemax] of the actual number of coefficients averaged at octave j. (these
%                       numbers are not powers of two because of the border effects at each scale). 
%             wtype:    The kind of weighting used
%                       0 -  no weigthing  (ie uniform weights)
%                       1 -  1/nj weights  (suitable for fully Gaussian data)
%                       2 -  use variance estimates varj
%             j1:       the lower limit of the scales chosen,  1<= j1 <= scalemax-1  
%             j2:       the upper limit of the octaves chosen, 2<= j2 <=
%             scalemax
%  Output:    slope:       estimate of the linear regression(s)
%             intercept:   intercept of the linear regression(s)               
%
% Herwig Wendt, Lyon, 2006 - 2008

jj = j1:j2 ;
J  = length(jj); 
njj  = nj(jj);
yjj  = yj(jj);
varyjj= varyj(jj);

if     wtype==0  % uniform weights
   wvarjj = ones(1,J);   
   wstr = 'Uniform';
elseif wtype==1  % Gaussian type weights
   wvarjj = 1./njj;        
   wstr = 'Gaussian';
elseif  wtype==2   % weights from data 
   wvarjj = varyjj; 
   wstr = 'Estimated';
else % all other cases
   fprintf('** Weight option not recognised, using uniform weights\n')
   wvarjj = ones(1,J); 
   wstr = 'Uniform';
end

% use weighted regression formula in all cases
S0 = sum(1./wvarjj) ;
S1 = sum(jj./wvarjj) ;
S2 = sum(jj.^2./wvarjj) ;
wjj = (S0 * jj - S1) ./ wvarjj / (S0*S2-S1*S1) ;
vjj = (S2 - S1 * jj) ./ wvarjj / (S0*S2-S1*S1) ;


%  Estimate  zeta
slope  = sum(wjj .* yjj ) ;       % zeta is just the slope, unbiased regardless of the weights
intercept     = sum(vjj .* yjj ) ;       % intercept  'a'
