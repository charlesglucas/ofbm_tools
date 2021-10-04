function [zetaest,aest]=Regrmat(yj,varyj,nj,wtype, j1,j2)
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
%  Output:    zetaest:   estimate of the linear regression(s)
%               aeta:    intercept of the linear regression(s)
%
% Herwig Wendt, Lyon, 2006 - 2008


%--- Initialise
[sz] = size(yj);
DIM=length(sz);
if (DIM==2)&(min(sz)==1)
    VECTOR=1;    
    if sz(1)>sz(2); yj=yj';[sz] = size(yj); end
else
    VECTOR=0;
end

scalemax=sz(end);
jr=1:scalemax;

%--- Check if only 1 estimate
if wtype~=2
    varyj=ones(size(yj));
end
if wtype~=1
    nj=1./jr; nj=nj/nj(end);
end

%--- Check and clean up the inputted j1 and j2
j1 = max(1,j1);                    % make sure j1 is not too small
j2 = max(j1+1,min(j2,scalemax));   % make sure j2 is not chosen too large

%--- Convenience variables in the chosen scaling range
jj = j1:j2 ;
Jr  = length(jj); 

% ---- MULTIPLE ESTIMATES
    njj  = squeeze(repmat(nj(jj),[1 1 sz(1:end-1)])); % 1st dimension: jr
    njj = permute(njj, [2:length(sz) 1]);               % jr is last dimension; size(njj)=size(yj)

    % select scales
    YJ=shiftdim(yj, length(size(yj))-1);
    SYJ=size(YJ);
    YJJ=YJ(jj,:);
    YJJ=reshape(YJJ, [Jr SYJ(2:end)]);
    yjj=shiftdim(YJJ,1);

    JJr=squeeze(repmat(jj,[1 1 sz(1:end-1)])); % 1st dimension: jr
    JJr = permute(JJr, [2:length(sz) 1]);       % jr is last dimension; size(JJr)=size(yj)

    % select scales
    VARYJ=shiftdim(varyj, length(size(varyj))-1);
    SVARYJ=size(VARYJ);
    VARYJJ=VARYJ(jj,:);
    VARYJJ=reshape(VARYJJ, [Jr SVARYJ(2:end)]);
    varyjj=shiftdim(VARYJJ,1);


    %%%% Regression,  VERSION DETERMINISTE
    %% define the effective varjj's used, depending on weight type
    if     wtype==0  % uniform weights
        wvarjj = ones(size(njj));
        wstr = 'w Unif';
    elseif wtype==1  % Gaussian type weights
        wvarjj = 1./njj;
        wstr = 'w Gauss';
    elseif  wtype==2   % weights from data
        wvarjj = varyjj;
        wstr = 'w Est';
    else % all other cases
        fprintf('** Weight option not recognised, using uniform weights\n')
        wvarjj = ones(size(njj));
        wstr = 'w Unif';
    end

    % use weighted regression formula in all cases
    S0 = sum(1./wvarjj,DIM) ;
    S1 = sum(JJr./wvarjj,DIM) ;
    S2 = sum(JJr.^2./wvarjj,DIM) ;
    S0=squeeze(repmat(S0,[ones(1, DIM) Jr]));
    S1=squeeze(repmat(S1,[ones(1, DIM) Jr]));
    S2=squeeze(repmat(S2,[ones(1, DIM) Jr]));

    wjj = (S0 .* JJr - S1) ./ wvarjj ./ (S0.*S2-S1.*S1) ;
    vjj = (S2 - S1 .* JJr) ./ wvarjj ./ (S0.*S2-S1.*S1) ;

    %  Estimate  zeta
    zetaest  = (sum(wjj .* yjj,DIM )) ;       % zeta is just the slope, unbiased regardless of the weights
    aest     = (sum(vjj .* yjj,DIM )) ;       % intercept  'a'

    % FORMAT OUTPUT
    zetaest=shiftdim(zetaest,1);
    aest=shiftdim(aest,1);

end