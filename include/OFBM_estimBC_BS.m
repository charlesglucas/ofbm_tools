
function [est,estbc] = OFBM_estimBC_BS(data,params) 

% returns the univariate-like and multivariate estimations in a structure
% 'est' and the bias corrected estimation in a structure 'estbc', which is
% empty when Jref = 0
%
%
%  Input:     data:   matrix of size P x N, with P number of components and N 
%                   the sample size, corresponding to an ofBm model     
%             params: structure of parameters
%                   - wtype: the kind of weighting used in the linear
%                   regressions:
%                       0 -  no weigthing  (ie uniform weights)
%                       1 -  1/nj weights  (suitable for fully Gaussian data)
%                       2 -  use variance of the estimates
%                   - Nwt: order of the Debauchies wavelet
%                   - Jref: reference scale under which several wavelet
%                   spectra are computed with the same number wavelet coeficients
%                   - j2: last scale for analysis
%                   - j1: first scale for analysis
%                   - NB: number of bootstrap resampling
%                   - LB: number of blocks for the bootstrap resampling
%                   - FigNum: number of the figures to plot, nothing is
%                   plotted if FigNum < 1
%                   - nbcoeffmin: coefficent to adjust the maximum scale of 
%                   analysis according to the sample sizes at each scale
%                   - nbcompmaxplot: the maximum number of components for
%                   which wavelet spectrum is plotted
%                   - Jrefplot: scale of shifting to overlap the structure
%                   function plots
%
%  Output:    est:     structure related to the univariate-like and classical
%                   multivariate estimation
%                   - JMj: P vectors indicating the maximum octaves for
%                   which wavelet spectrum rank is not defficient for each
%                   component
%                   - j2j: P vectors indicating the maximum octave of
%                   analysis accounting for the rank defficiency for each
%                   component
%                   - Nj: NJ vector indicating the number of available wavelet
%                   coefficients at each scale
%                   - WW: P x P x NJ wavelet spectrum across the NJ
%                   scales such that P < N   
%                   - CondNumber: NJ condition numbers of wavelet spectrum
%                   across scales
%                   - rankS: NJ ranks of wavelet spectrum across scales
%                   - wavcov: P x P x NJ wavelet coherence with NJ the number of
%                   scales such that P < N  
%                   - lambdaj: NJ x P matrix containing the P eigenvalues of the
%                   wavelet spectrum WW at each scale
%                   - VVj: P x P x NJ matrix containing the P eigenvectors of the
%                   wavelet spectrum WW at each scale
%                   - lambda : P multivariate estimates of Hurst exponents,
%                   up to an   additive constant
%                   - lambdasort : P sorted multivariate estimates of Hurst exponents,
%                   up to an additive constant
%                   - h: P multivariate estimates of Hurst exponents
%                   - lambdaBS : NB x P matrix containing the NB bootstrap resamples of 
%                   multivariate estimates of Hurst exponents, up to an additive constant
%                   - lambdajBS: NJ x NB x P matrix containing the NB bootstrap resambles
%                   of the P eigenvalues of the wavelet spectrum WW at each scale
%                   - hBS: NB x P matrix containing the NB bootstrap resamples of multivariate 
%                   Hurst exponent estimates
%                   - intercept: P intercepts from the linear regressions of
%                   wavelet spectrum log-eigenvalues against scales
%                   - hU: P x P univariate estimates and cross-estimates of Hurst exponents
%                   - interceptU: P intercepts from the linear regressions of
%                   wavelet spectrum diagonal coefficients against scales
%             estbc:  structure related to a bias corrected multivarate estimation,
%                   which is empty when Jref = 0
%                   - WW: NJ cells containing P x P x Nj wavelet spectrum
%                   across scales such that P < N
%                   - lambdaj: NJ x P matrix containing the P corrected eigenvalues of the
%                   wavelet spectrum WW at each scale
%                   - lambda : P bias corrected multivariate estimates of Hurst exponents,
%                   up to an additive constant
%                   - lambdasort : P sorted bias corrected multivariate estimates of Hurst 
%                   exponents, up to an additive constant
%                   - h: P bias corrected multivariate estimates of Hurst exponents
%                   - lambdajBS: NJ x NB x P matrix containing NB resamples of the P corrected
%                   eigenvalues of the wavelet spectrum WW at each scale
%                   - lambdaBS : NB x P matrix containing the NB bootstrap resamples of 
%                   bias corrected multivariate estimates of Hurst exponents, up to an 
%                   additive constant
%                   - hBS: NB x P matrix containing the NB bootstrap resamples of bias corrected 
%                   multivariate Hurst exponent estimates
%                   - intercept: P intercepts from the linear regressions of corrected wavelet
%                   spectrum log-eigenvalues against scales
%
% Herwig Wendt, Patrice Abry and Charles-Gérard Lucas, ENS Lyon, 2020 - 2022


estbc = {};

if ~isfield(params,'FBM'), params.FBM=1; end
if ~isfield(params,'wtype'), params.wtype=1; end
if ~isfield(params,'Nwt'), params.Nwt=2; end
if ~isfield(params,'Jref'), params.Jref=params.j2; end
if ~isfield(params,'j2'), params.j2=log2(size(data,2)) - params.Nwt - 3 ; end
if ~isfield(params,'j1'), params.j1=params.j2-2; end
if ~isfield(params,'NB'), params.NB=0; end
if ~isfield(params,'LB'), params.LB=1*params.Nwt; end
if ~isfield(params,'FigNum'), params.FigNum=0; end
if ~isfield(params,'nbcoeffmin'), params.nbcoeffmin = 8 ; end
if ~isfield(params,'nbcompmaxplot'), params.nbcompmaxplot = 5 ; end
if ~isfield(params,'Jrefplot'), params.Jrefplot = 1 ; end
if ~isfield(params,'gamint'), params.gamint = 0.5; end
if params.j2 <= params.j1, params.j1 = params.j2-2; end

Nwt = params.Nwt ; 
j1 = params.j1 ;
j2 = params.j2 ; 
wtype = params.wtype ;
Jref = params.Jref ; 
NB = params.NB ; 
LB = params.LB ; 
FigNum = params.FigNum ;
nbcoeffmin = params.nbcoeffmin ; 
Jrefplot = params.Jrefplot ; 
gamint = params.gamint;

markersize = 15 ; 
markersize2 = 6 ; 
linewidth = 2 ; 
linewidth2 = 1 ; 
fontsize = 20 ;
fontsize2 = 12 ;

[P,nn] = size(data) ; 
if P>nn
    tmp = nn ; 
    nn = P ;
    P = tmp ; 
end

for k = 1:1:P
    [cx,lx, nj] = DxLx1d(data(k,:),Nwt,gamint);
     NJ=nj.W;
     j2 = min(j2,length(NJ)) ; 
     WD{k} = cx ; 
end

JM = length(nj.W) ; JM=JM(end) ;
JMtmp=find(nj.W>=P) ; JMtmp=JMtmp(end)  ;

IID = 0;
% --- bootstrap
if NB
    WD=Dx_bootstrap_OFBM(WD,NB,LB,j1,IID); 
end
J1BS=j1;
  
JM=find(nj.W>=P); JM=JM(end); j2=min(j2,JM); 
nj.W=nj.W(1:JM);
nj.L=nj.L(1:JM);
NJ=nj.W;

for p=P:-1:1
    JMtmp2=find(nj.W>=P-p+nbcoeffmin) ; 
    JMj(p) = JMtmp2(end)  ; 
    j2j(p) = min(j2,JMj(p)) ; 
end
est.JMj = JMj ; 
est.j2j = j2j ; 

Jref = min([Jref,j2j]);

% if Jref==0, don't divide and use standard estimation
if Jref ~=0
    Ldiv = nj.W(Jref);
    NJ2 = max(Ldiv,floor(nj.W/Ldiv)*Ldiv);
    NJ2 = min([NJ2;nj.W]);
end


for k = 1:1:P
    for m = 1:1:P
        for j=1:length(NJ)
            tmp=mean(real(WD{k}(j).value_noabs.*WD{m}(j).value_noabs));
            WW(k,m,j) = tmp ; 
            % --- bootstrap
            if j>=J1BS && NB
                tmpBS=mean(real(WD{k}(j).value_noabs_bs.*WD{m}(j).value_noabs_bs),2)';
                WWBS(k,m,:,j) = tmpBS ;
            end
            
            if Jref ~= 0    
                nntmp=min(Ldiv,NJ2(j));
                tmp=mean(reshape(real(WD{k}(j).value_noabs(1:NJ2(j)).*WD{m}(j).value_noabs(1:NJ2(j))),nntmp,[]));
                WWbc{j}(k,m,:) = tmp ; 
                % --- bootstrap
                if j>=J1BS && NB
                    tmpBS=mean(reshape(real(WD{k}(j).value_noabs_bs(:,1:NJ2(j)).*WD{m}(j).value_noabs_bs(:,1:NJ2(j))), NB, nntmp,[]),2 );
                    WWBSbc{j}(k,m,:,:) = tmpBS ;
                end
            end
        end
    end
end


[a,b,J] = size(WW) ; 
if Jref ~= 0, Jbc = length(WWbc) ; end

for j = 1:JM
    tmp2 = WW(:,:,j) ;
    CondNumber(j) = rcond(tmp2) ; 
    rankS(j) = rank(tmp2) ; 
    [Vtmp,Dtmp] = eig(tmp2) ;
    for k = 1:1:P
        lambda(j,k) = Dtmp(k,k) ;
    end
    VV(:,:,j) = Vtmp ;
    
    if Jref ~= 0
        % trim size of S and sum
        Vtmpbc=[]; lamb=[];
        tmp1 = WWbc{j}; [~,~,nb]=size(tmp1);
        clear lamb Vtmp
        for ib=1:nb
            tmp2=squeeze(WWbc{j}(:,:,ib));
            [V0,Dtmpbc]=eig(tmp2);
            lamb(ib,:)=diag(Dtmpbc);
            lamb(ib,lamb(ib,:)<=0)=NaN;
            Vtmpbc(ib,:,:)=V0;
        end
        if nb>1
            lambdabc(j,:) = 2.^nanmean(log2(lamb));
            %lambdabc(j,:) = nanmean(lamb);
            %lambdabc(j,:) = 2.^nanmedian(log2(lamb));
            VVbc(:,:,j) = squeeze(mean(Vtmpbc)) ;
        else
            lambdabc(j,:) = lamb;
            VVbc(:,:,j) = Vtmpbc ;
        end
    end
end

% --- bootstrap
for nb=1:NB
    for j = J1BS:JM
        tmp = WWBS(:,:,nb,j) ;
        [Vtmp,Dtmp] = eig(tmp) ;
        for k = 1:1:P
            lambdaBS(k,nb,j) = Dtmp(k,k) ;
        end
    end
end

if Jref ~= 0
    for nb=1:NB
        for j = J1BS:JM
            tmp1 = WWBSbc{j}; [~,~,~,ndiv]=size(tmp1);
            clear lambBSbc
            for ib=1:ndiv
                tmp = WWBSbc{j}(:,:,nb,ib) ;
                [Vtmpbc,Dtmpbc] = eig(tmp) ;
                lambBSbc(:,ib)=diag(Dtmpbc);
                lambBSbc(lambBSbc(:,ib)<=0,ib)=NaN;
            end
            try
            if ndiv>1
                lambdaBSbc(:,nb,j) = 2.^squeeze(mean(log2(lambBSbc),2));
            else
                lambdaBSbc(:,nb,j) = lambBSbc;
            end
            catch
                keyboard
            end
        end
    end
end

for k = 1:1:P
    for m = 1:1:P
        wavcov(k,m,:)= squeeze(WW(k,m,:))./sqrt(squeeze(WW(k,k,:)).*squeeze(WW(m,m,:))) ;
    end
end

for nb=1:NB
    for k = 1:1:P
        for m = 1:1:P
            wavcovBS(k,m,nb,:)= squeeze(WWBS(k,m,nb,:))./sqrt(squeeze(WWBS(k,k,nb,:)).*squeeze(WWBS(m,m,nb,:))) ;
        end
    end
end

if Jref ~= 0
    for j = 1:JM
        tmp1 = WWbc{j}; [~,~,nb]=size(tmp1);
        wavcovjb=[];
        for k = 1:1:P
            for m = 1:1:P
                wavcovjb(:,k,m,j)=squeeze(WWbc{j}(k,m,:))./sqrt(squeeze(WWbc{j}(k,k,:)).*squeeze(WWbc{j}(m,m,:))) ;
                wavcovbc(k,m,j)=squeeze(mean(wavcovjb(:,k,m,j)));
            end
        end
        
        meanWWbc(:,:,j) = mean(WWbc{j},3);
    end
end

JJ = 1:1:J ; 
JJJ = [j1:1:j2] ; 
est.lambdaj = lambda ; 

if Jref ~= 0
    JJbc = 1:1:Jbc ; 
    estbc.lambdaj = lambdabc;
    
end

est.VVj = VV ;
est.wavcov = wavcov ; 
est.VV = squeeze(VV(:,:,J-2)) ; 
est.beta = est.VV(1,2)/est.VV(2,2) ; 
est.WW = WW ;
for j=1:1:J
    est.betaj(j) = VV(1,2,j)/VV(2,2,j) ;
end
est.beta2 = mean(est.betaj(J-1:J)) ;   
est.Nj = NJ ; 
est.CondNumber = CondNumber ; 
est.RankS = rankS ;

if NB, est.lambdajBS = lambdaBS; end
[est.hU,est.interceptU]=Regrmat(log2(squeeze(abs(est.WW))),ones(P,P,length(NJ)),NJ,wtype, j1,j2);
est.hU = est.hU/2;


if Jref ~= 0
    estbc.WW = WWbc ;
    if NB, estbc.lambdajBS = lambdaBSbc; end
end

for p=1:1:P
    tmp = log2(lambda(:,p)') ; 
    [slope,intercept] = SmartLDplot(tmp(1:j2j(p)),NJ(1:j2j(p)),ones(1,j2j(p)),j1,j2j(p),wtype) ;
    est.lambda(p) = (slope)/2 ;  
    est.intercept(p) = intercept;
    
    % --- bootstrap
    if NB
        tmpBS=log2(abs(squeeze(lambdaBS(p,:,:))));
        [slopeBS,Vzeta,Q,interceptBS]=MFA_BS_regrmat(tmpBS(:,1:j2j(p)),ones(size(tmpBS,1),j2j(p)),NJ(1:j2j(p)),wtype, j1,j2j(p));
        est.lambdaBS(:,p) = (slopeBS)/2 ;  
        est.hBS(:,p) = (slopeBS)/2  + 1/2 - params.FBM ;  
    end

end
est.h = est.lambda + 1/2 - params.FBM ;

[est.lambdasort,~]=sort(est.lambda,'ascend');

if Jref ~= 0
    for p=1:1:P
        tmp = log2(lambdabc(:,p)') ; 
        [slope,intercept] = SmartLDplot(tmp(1:Jref),NJ(1:Jref),ones(1,Jref),j1,Jref,wtype) ; 
        estbc.lambda(p) = (slope)/2 ;  
        estbc.intercept(p) = intercept;

        % --- bootstrap
        if NB
            tmpBS=log2(abs(squeeze(lambdaBSbc(p,:,:))));
            [slopeBS,Vzeta,Q,interceptBS]=MFA_BS_regrmat(tmpBS(:,1:Jref),ones(size(tmpBS,1),Jref),NJ(1:Jref),wtype,j1,Jref);
            estbc.lambdaBS(:,p) = (slopeBS)/2 ;  
            estbc.hBS(:,p) = (slopeBS)/2 + 1/2 - params.FBM ;  
        end
    end

    [estbc.lambdasort,~]=sort(estbc.lambda,'ascend');
    estbc.h = estbc.lambda + 1/2 - params.FBM ;
end

if FigNum > 0 
    JMM = length(lambda(:,P)) ; 
    if P<=params.nbcompmaxplot
        figure(FigNum) ; clf 
        for p = 1:1:P
            for m = 1:1:P
                if NB
                    ic = 1.96/sqrt(NB)*std(log2(abs(squeeze(WWBS))),[],3);
                    ic2 = 1.96/sqrt(NB)*std(log2(abs(squeeze(wavcovBS))),[],3);
                end
                subplot(P,P,(p-1)*P+m); 
                if p~=m
                    plot(squeeze(est.wavcov(p,m,:)),'ob-','MarkerSize',markersize2,'LineWidth',linewidth2) ;
                    if NB, errorbar(squeeze(est.wavcov(p,m,:)),squeeze(ic2(p,m,:)),'ob-','MarkerSize',markersize2,'LineWidth',linewidth2) ; end
                    grid on ; 
                    axis([0 JMM+1 -1.01 1.01]) ; hold on ;
                    plot(zeros(size(squeeze(est.wavcov(p,m,:)))),'k--','MarkerSize',markersize2,'LineWidth',linewidth2) 
                    xlabel('$j = \log_2 2^j$','Interpreter','Latex','FontSize',fontsize2) ;
                    ylabel(['$C_{',num2str(p),num2str(m),'}(2^j)$'],'Interpreter','Latex','FontSize',fontsize2)
                    pbaspect([1.1,1,1])
                    set(gca,'LineWidth',linewidth2)
                else
                    plot(log2(abs(squeeze(est.WW(p,m,:)))),'ob-','MarkerSize',markersize2,'LineWidth',linewidth2) ; 
                    if NB, errorbar(log2(abs(squeeze(est.WW(p,m,:)))),squeeze(ic(p,m,:)),'ob-','MarkerSize',markersize2,'LineWidth',linewidth2); end
                    grid on ; 
                    V((p-1)*P+m,:) = axis ; 
                    xlabel('$j = \log_2 2^j$','Interpreter','Latex','FontSize',fontsize2) ;
                    ylabel(['$\log_2 S_{',num2str(p),num2str(p),'}(2^j)$'],'Interpreter','Latex','FontSize',fontsize2)
                    pbaspect([1.1,1,1])
                    set(gca,'LineWidth',linewidth2)
                end
            end
        end
        clear Vf
        Vf(1) = 0 ; Vf(2) = JMM+1 ; Vf(3) = min(V(:,3)) ; Vf(4) = max(V(:,4)) ; 
        for p = 1:1:P, subplot(P,P,(p-1)*P+p) ; axis(Vf)  ; end
        fig.Position = [57 160 834 643];
        sgtitle('Wavelet spectrum','Interpreter','Latex')
    end

    for j = 1:JM, dWW(j,:) = diag(WW(:,:,j)); end
    for nb=1:NB, for j = 1:JM, dWWBS(j,:,nb) = diag(WWBS(:,:,nb,j)); end; end
    if NB, ic = 1.96/sqrt(NB)*std(log2(abs(dWWBS)),[],3); end
    hh = figure(FigNum+1) ; clf ; 
    set(gca,'FontSize',fontsize) ;   
    for p = 1:1:P
        JJ =1:JMM ;
        JJJ = j1:j2 ;
        plot(JJ,log2(abs(dWW(:,p))),'ok-','MarkerSize',markersize,'LineWidth',linewidth) ;
        if NB, errorbar(JJ,log2(abs(dWW(:,p))),ic(:,p),'ok-','MarkerSize',markersize,'LineWidth',linewidth) ; end
        grid on ; hold on ;
        pp = polyfit(JJJ,log2(abs(dWW(JJJ,p)))',1) ; 
        qq = polyval(pp,JJ) ; 
        plot(1:JMM,qq,'k--','MarkerSize',markersize,'LineWidth',linewidth) ; 
        plot(JJJ,qq(JJJ),'k-','MarkerSize',markersize,'LineWidth',linewidth) ; 
        Vax = axis ; plot([JMtmp JMtmp],[Vax(3) Vax(4)],'k--') ; 
    end
    xlabel('$j = \log_2 2^j$','Interpreter','Latex','FontSize',fontsize) ;
    ylabel('$\log_2 S_{mm}(2^j)$','Interpreter','Latex','FontSize',fontsize)
    set(gca,'FontSize',fontsize,'TickLabelInterpreter','Latex')

    if NB, ic = 1.96/sqrt(NB)*squeeze(std(log2(lambdaBS),[],2))'; end
    hh = figure(FigNum+2) ; clf ; 
    set(gca,'FontSize',fontsize) ; 
    for p = 1:1:P
        JJ =1:JMj(p) ;
        JJJ = j1:j2j(p) ;
        plot(JJ,log2(lambda(JJ,p)),'ob-','MarkerSize',markersize,'LineWidth',linewidth) ;
        if NB, errorbar(JJ,log2(lambda(JJ,p)),ic(JJ,p),'ob-','MarkerSize',markersize,'LineWidth',linewidth) ; end 
        grid on ; hold on ;
        pp = polyfit(JJJ,log2(abs(lambda(JJJ,p)))',1) ; 
        qq = polyval(pp,1:JMM) ; 
        plot(1:JMM,qq,'k--','MarkerSize',markersize,'LineWidth',linewidth) ; 
        plot(JJJ,qq(JJJ),'k-','MarkerSize',markersize,'LineWidth',linewidth) ; 
        Vax = axis ; plot([JMtmp JMtmp],[Vax(3) Vax(4)],'k--') ; 
    end
    xlabel('$j = \log_2 2^j$','Interpreter','Latex','FontSize',fontsize) ;
    ylabel('$\log_2 \lambda_m(2^j)$','Interpreter','Latex','FontSize',fontsize)
    set(gca,'FontSize',fontsize,'TickLabelInterpreter','Latex')

    if Jref==0
        hh = figure(FigNum+3) ; clf ; 
        set(gca,'FontSize',fontsize) ; 
        for p = 1:1:P
            JJ =1:JMj(p) ;
            plot(JJ,log2(lambda(JJ,p)),'ob-','MarkerSize',markersize,'LineWidth',linewidth) ; grid on ; hold on 
            plot(log2(abs(squeeze(est.WW(p,p,JJ)))),'ok-','MarkerSize',markersize,'LineWidth',linewidth) ; 
        end
        xlabel('$j = \log_2 2^j$','Interpreter','Latex','FontSize',fontsize) ;
        ylabel('$\log_2 \lambda_m(2^j)$','Interpreter','Latex','FontSize',fontsize)
        set(gca,'FontSize',fontsize,'TickLabelInterpreter','Latex')
    else
        if NB, ic = 1.96/sqrt(NB)*squeeze(std(log2(lambdaBSbc),[],2))'; end
        hh = figure(FigNum+3) ; clf ; 
        set(gca,'FontSize',fontsize) ; 
        for p = 1:1:P
            JJ =1:JMj(p) ;
            JJJ = j1:j2j(p) ;
            plot(JJ,log2(lambdabc(JJ,p)),'or-','MarkerSize',markersize,'LineWidth',linewidth) ; 
            if NB, errorbar(JJ,log2(lambdabc(JJ,p)),ic(JJ,p),'or-','MarkerSize',markersize,'LineWidth',linewidth) ; end
            grid on ; hold on ;
            pp = polyfit(JJJ,log2(abs(lambdabc(JJJ,p)))',1) ; 
            qq = polyval(pp,1:JMM) ; 
            plot(1:JMM,qq,'k--','MarkerSize',markersize,'LineWidth',linewidth) ; 
            plot(JJJ,qq(JJJ),'k-','MarkerSize',markersize,'LineWidth',linewidth) ; 
            Vax = axis ; plot([JMtmp JMtmp],[Vax(3) Vax(4)],'k--') ; 
        end
        xlabel('$j = \log_2 2^j$','Interpreter','Latex','FontSize',fontsize) ;
        ylabel('$\log_2 \lambda_m(2^j)$','Interpreter','Latex','FontSize',fontsize)
        set(gca,'FontSize',fontsize,'TickLabelInterpreter','Latex')

        hh = figure(FigNum+4) ; clf ; 
        set(gca,'FontSize',fontsize) ; 
        for p = 1:1:P
            JJ =1:JMj(p) ;
            plot(log2(abs(squeeze(est.WW(p,p,JJ)))),'ok-','MarkerSize',markersize,'LineWidth',linewidth) ;  grid on ; hold on ;
            plot(JJ,log2(abs(lambda(JJ,p))),'ob-','MarkerSize',markersize,'LineWidth',linewidth) ;
            plot(JJ,log2(abs(lambdabc(JJ,p))),'or-','MarkerSize',markersize,'LineWidth',linewidth) ;
        end
        xlabel('$j = \log_2 2^j$','Interpreter','Latex','FontSize',fontsize) ;
        ylabel('$\log_2 \lambda_m(2^j)$','Interpreter','Latex','FontSize',fontsize)
        set(gca,'FontSize',fontsize,'TickLabelInterpreter','Latex')

        hh = figure(FigNum+5) ; clf ; 
        set(gca,'FontSize',fontsize) ; 
        for p = 1:1:P
            JJ =1:JMj(p) ;
            plot(log2(abs(squeeze(est.WW(p,p,JJ))))-log2(abs(squeeze(est.WW(p,p,Jrefplot)))),'ok-','MarkerSize',markersize,'LineWidth',linewidth) ;  grid on ; hold on ;
            plot(JJ,log2(abs(lambda(JJ,p)))-log2(abs(lambda(Jrefplot,p))),'ob-','MarkerSize',markersize,'LineWidth',linewidth) ;
            plot(JJ,log2(abs(lambdabc(JJ,p)))-log2(abs(lambdabc(Jrefplot,p))),'or-','MarkerSize',markersize,'LineWidth',linewidth) ;
        end
        xlabel('$j = \log_2 2^j$','Interpreter','Latex','FontSize',fontsize) ;
        ylabel('$\log_2 \lambda_m(2^j)  - \log_2 \lambda_m(2^{\textrm{Jrefplot}})$','Interpreter','Latex','FontSize',fontsize)
        set(gca,'FontSize',fontsize,'TickLabelInterpreter','Latex')
    end
end

end