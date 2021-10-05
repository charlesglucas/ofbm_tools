

% function [est,estbc,nj,JM] = OFBM_estim_divBC_BS(data,Nwt,j1,j2,wtype,Jref,NB,LB,FigNum)
function [est,estbc] = OFBM_estimBC_BS(data,params) 
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
%                   which wavalet spectrum is plotted
%                   - Jrefplot: scale of shifting to overlap the structure
%                   function plots
%
%  Output:    est:     structure related to the univariate-like and classical
%                   multivariate estimation
%             estbc:   structure related to a bias corrected multivarate estimation,
%                   which is empty when Jref = 0
%

estbc = {};

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


if params.Jref > params.j2, params.Jref = params.j2; end
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

gamint = 0.5; 
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


Jref=min(Jref,j2);


% HW: if Jref==0, don't divide and use standard estimation
if Jref ~=0
    Ldiv = nj.W(Jref);
    NJ2=max(Ldiv,floor(nj.W/Ldiv)*Ldiv);
    NJ2=min([NJ2;nj.W]);
end

for k = 1:1:P
    for m = 1:1:P
        for j=1:length(NJ)
            tmp=mean(real(WD{k}(j).value_noabs.*WD{m}(j).value_noabs));
            WW(k,m,j) = tmp ; 
            % --- bootstrap
            if j>=J1BS & NB
                tmpBS=mean(real(WD{k}(j).value_noabs_bs.*WD{m}(j).value_noabs_bs),2)';
                WWBS(k,m,:,j) = tmpBS ;
            end
            
            if Jref ~= 0
                nntmp=min(Ldiv,NJ2(j));
                tmp=mean(reshape(real(WD{k}(j).value_noabs(1:NJ2(j)).*WD{m}(j).value_noabs(1:NJ2(j))),nntmp,[]));
                WWbc{j}(k,m,:) = tmp ; 
                % --- bootstrap
                if j>=J1BS & NB
                    tmpBS=mean(reshape( real(WD{k}(j).value_noabs_bs(:,1:NJ2(j)).*WD{m}(j).value_noabs_bs(:,1:NJ2(j))), NB, nntmp,[]),2 );
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
            [V0,Dtmpbc] = eig(tmp2) ;
            lamb(ib,:)=diag(Dtmpbc);
            Vtmpbc(ib,:,:)=V0;
        end
        if nb>1
            lambdabc(j,:)=mean(lamb);
            VVbc(:,:,j) = squeeze(mean(Vtmpbc)) ;
        else
            lambdabc(j,:)=lamb;
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
%            VVBS(:,:,nb,j) = Vtmp ;
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
            end
            try
            if ndiv>1
                lambdaBSbc(:,nb,j) = squeeze(mean(lambBSbc,2));
            else
                lambdaBSbc(:,nb,j) = lambBSbc;
            end
            catch
                keyboard
            end
            %VVBS(:,:,nb,j) = Vtmp ;
        end
    end
end

% --- bootstrap
% for nb=1:NB
%     for j = J1BS:JM
%         tmp = WWBS(:,:,nb,j) ;
%         [Vtmp,Dtmp] = eig(tmp) ;
%         for k = 1:1:P
%             lambdaBS(k,nb,j) = Dtmp(k,k) ;
%         end
%         VVBS(:,:,nb,j) = Vtmp ;
%     end
% end
% if Bcorr==1;
%     for j = J1BS:JM
%         for k = 1:1:P
%             lambda(j,k) = 2*lambda(j,k) - mean(squeeze(lambdaBS(k,:,j)));
%         end
%     end
% end

%keyboard
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
    estbc.lambdaj = lambdabc ; 
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
est.nj = nj ; 
est.JM = JMtmp ; 
est.CondNumber = CondNumber ; 
est.RankS = rankS ;
est.j1 = j1 ; 
est.j2 = j2 ;
est.wtype = wtype ;

% --- bootstrap
% est.lambdajmbs = mean(squeeze(lambdaBS(k,:,j))) ; 
% est.WWBS = WWBS ;
if NB, est.lambdajBS = lambdaBS; end
est.hU=Regrmat(log2(squeeze(abs(est.WW))),ones(P,P,length(NJ)),NJ,wtype, j1,j2)/2 + 1/2 - params.FBM;


if Jref ~= 0
    estbc.VVj = VVbc ;
    estbc.wavcov = wavcovbc ; 
    estbc.VV = squeeze(VVbc(:,:,Jbc-2)) ; 
    estbc.beta = estbc.VV(1,2)/estbc.VV(2,2) ; 
    estbc.WW = WWbc ;
    for j=1:1:Jbc
        estbc.betaj(j) = VVbc(1,2,j)/VVbc(2,2,j) ;
    end
    estbc.beta2 = mean(estbc.betaj(Jbc-1:Jbc)) ;  
    if NB, estbc.lambdajBS = lambdaBSbc; end
end

for p=1:1:P
    tmp = log2(lambda(:,p)') ; 
     [slope,intercept] = SmartLDplot(tmp,NJ,ones(size(tmp)),j1,j2,wtype) ; 
    est.lambda(p) = (slope)/2 ;  
    est.intercept(p) = intercept;
    
    % --- bootstrap
    if NB
        tmpBS=log2(abs(squeeze(lambdaBS(p,:,:))));
        [slopeBS,Vzeta,Q,interceptBS]=MFA_BS_regrmat(tmpBS,ones(size(tmpBS)),NJ,wtype, j1,j2);
        est.lambdaBS(:,p) = (slopeBS)/2 ;  
        est.meanBS(p)=mean(est.lambdaBS(:,p));
        est.varBS(p)=var(est.lambdaBS(:,p));
    end

end
est.h = est.lambda + 1/2 - params.FBM ;
[~,est.sortid]=sort(est.lambda,'ascend');

[est.lambdasort,sortid]=sort(est.lambda,'ascend');

if Jref ~= 0
    for p=1:1:P
        tmp = log2(lambdabc(:,p)') ; 
         [slope,intercept] = SmartLDplot(tmp,NJ,ones(size(tmp)),j1,j2,wtype) ; 
        estbc.lambda(p) = (slope)/2 ;  
        estbc.intercept(p) = intercept;

        % --- bootstrap
        if NB
            tmpBS=log2(abs(squeeze(lambdaBSbc(p,:,:))));
            [slopeBS,Vzeta,Q,interceptBS]=MFA_BS_regrmat(tmpBS,ones(size(tmpBS)),NJ,wtype, j1,j2);
            estbc.lambdaBS(:,p) = (slopeBS)/2 ;  
            estbc.meanBS(p)=mean(estbc.lambdaBS(:,p));
            estbc.varBS(p)=var(estbc.lambdaBS(:,p));
        end
    end

    [estbc.lambdasort,sortid]=sort(estbc.lambda,'ascend');
    estbc.h = estbc.lambda + 1/2 - params.FBM ;
    [~,estbc.sortid]=sort(estbc.lambda,'ascend');
end

% --- bootstrap
% est.lambdaBSsort=sort(est.lambdaBS,2,'ascend');
% for p1=1:1:P
%     %est.lambda(p1)=2*est.lambda(p1)-mean(est.lambdaBS(:,p1));
%     est.meanBS(p1)=mean(est.lambdaBS(:,p1));
%     est.varBS(p1)=var(est.lambdaBS(:,p1));
%     est.meanBSsort(p1)=mean(est.lambdaBSsort(:,p1));
%     est.varBSsort(p1)=var(est.lambdaBSsort(:,p1));
% end

if FigNum > 0 
   
    
    JMM = length(lambda(:,P)) ; 
    %t = tiledlayout(P,P);
    %t.TileSpacing = 'tight';
if params.P<=params.nbcompmaxplot
    figure(FigNum) ; clf ;
    for p = 1:1:P
        for m = 1:1:P
            if NB
                ic = 1.96/sqrt(NB)*std(log2(abs(squeeze(WWBS))),[],3);
                ic2 = 1.96/sqrt(NB)*std(log2(abs(squeeze(wavcovBS))),[],3);
            end
            subplot(P,P,(p-1)*P+m) ; 
            %nexttile
            if p~=m
                plot(squeeze(est.wavcov(p,m,:)),'ob-','MarkerSize',markersize2,'LineWidth',linewidth2) ;
                if NB, errorbar(squeeze(est.wavcov(p,m,:)),squeeze(ic2(p,m,:)),'ob-','MarkerSize',markersize2,'LineWidth',linewidth2) ; end
                grid on ; 
                axis([0 JMM+1 -1.01 1.01]) ; hold on ;
                plot(zeros(size(squeeze(est.wavcov(p,m,:)))),'k--','MarkerSize',markersize2,'LineWidth',linewidth2) 
                xlabel('$j = \log_2 2^j$','Interpreter','Latex','FontSize',fontsize2) ;
                ylabel(['$C_{',num2str(p),num2str(m),'}(2^j)$'],'Interpreter','Latex','FontSize',fontsize2)
                pbaspect([1.1,1,1])
            else
                plot(log2(abs(squeeze(est.WW(p,m,:)))),'ob-','MarkerSize',markersize2,'LineWidth',linewidth2) ; 
                if NB, errorbar(log2(abs(squeeze(est.WW(p,m,:)))),squeeze(ic(p,m,:)),'ob-','MarkerSize',markersize2,'LineWidth',linewidth2); end
                grid on ; 
                V((p-1)*P+m,:) = axis ; 
                xlabel('$j = \log_2 2^j$','Interpreter','Latex','FontSize',fontsize2) ;
                ylabel(['$S_{',num2str(p),num2str(p),'}(2^j)$'],'Interpreter','Latex','FontSize',fontsize2)
                pbaspect([1.1,1,1])
            end
        end
    end
    clear Vf
    Vf(1) = 0 ; Vf(2) = JMM+1 ; Vf(3) = min(V(:,3)) ; Vf(4) = max(V(:,4)) ; 
    for p = 1:1:P, subplot(P,P,(p-1)*P+p) ; axis(Vf)  ; end
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
    end
    
    hh = figure(FigNum+5) ; clf ; 
    set(gca,'FontSize',fontsize) ; 
    for p = 1:1:P
        JJ =1:JMj(p) ;
        plot(log2(abs(squeeze(est.WW(p,p,JJ))))-log2(abs(squeeze(est.WW(p,p,Jrefplot)))),'ok-','MarkerSize',markersize,'LineWidth',linewidth) ;  grid on ; hold on ;
        plot(JJ,log2(abs(lambda(JJ,p)))-log2(abs(lambda(Jrefplot,p))),'ob-','MarkerSize',markersize,'LineWidth',linewidth) ;
        plot(JJ,log2(abs(lambdabc(JJ,p)))-log2(abs(lambdabc(Jrefplot,p))),'or-','MarkerSize',markersize,'LineWidth',linewidth) ;
    end
    xlabel('$j = \log_2 2^j$','Interpreter','Latex','FontSize',fontsize) ;
    ylabel('$\log_2 \lambda_m(2^j) - \log_2 \lambda_m(2^{\textrm{Jrefplot}})$','Interpreter','Latex','FontSize',fontsize)
    set(gca,'FontSize',fontsize,'TickLabelInterpreter','Latex')
end
    
%    if Jref ~=0
%         JMM = length(lambdabc(:,P)) ; 
%         figure(FigNum+100) ; clf 
%         if P<11
%         for p = 1:1:P
%             for m = 1:1:P
%                 subplot(P,P,(p-1)*P+m) ; 
%                 if p~=m
%                     plot(squeeze(estbc.wavcov(p,m,:)),'or-') ; grid on ;  
%                     axis([0 JMM+1 -1.01 1.01]) ; hold on ;
%                     plot(zeros(size(squeeze(estbc.wavcov(p,m,:)))),'k--') 
%                 else
%                     plot(log2(abs(squeeze(meanWWbc(p,m,:)))),'or-') ; grid on ; 
%                     V((p-1)*P+m,:) = axis ; 
%                 end
%             end
%         end
%         clear Vf
%         Vf(1) = 0 ; Vf(2) = JMM+1 ; Vf(3) = min(V(:,3)) ; Vf(4) = max(V(:,4)) ; 
%         for p = 1:1:P, subplot(P,P,(p-1)*P+p) ; axis(Vf)  ; end
%         end
%    end

end