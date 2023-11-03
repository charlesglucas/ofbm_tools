%% Compare univariate multivariate self-similarity exponents
% CGL, October 2021.

clc
clear all
close all 
format compact

addpath('../include/')
load('../data/result_estimbc_sizeH6.mat')
P = size(data,1);

%% Run analysis
paramsEst.j1 = 5;
paramsEst.j2 = 10;
paramsEst.Jrefplot = 1;

[est,estbc] = OFBM_estimBC_BS(data,paramsEst);

%% Wavelet spectrum
markersize2 = 6 ; 
linewidth2 = 1 ;
fontsize2 = 12 ;

JMM = size(estbc.lambdaj,1) ; 
fig = figure(1) ; clf 
for p = 1:1:P
    for m = 1:1:P
        subplot(P,P,(p-1)*P+m); 
        if p~=m
            plot(squeeze(est.wavcov(p,m,:)),'ok-','MarkerSize',markersize2,'LineWidth',linewidth2) ;
            grid on ; 
            axis([0 JMM+1 -1.01 1.01]) ; hold on ;
            plot(zeros(size(squeeze(est.wavcov(p,m,:)))),'k--','MarkerSize',markersize2,'LineWidth',linewidth2) 
            xlabel('$j = \log_2 2^j$','Interpreter','Latex','FontSize',fontsize2) ;
            ylabel(['$C_{',num2str(p),num2str(m),'}(2^j)$'],'Interpreter','Latex','FontSize',fontsize2)
            pbaspect([1.1,1,1])
            set(gca,'LineWidth',linewidth2)
            xticks(1:4:JMM)
        else
            plot(log2(abs(squeeze(est.WW(p,m,:)))),'ob-','MarkerSize',markersize2,'LineWidth',linewidth2) ; 
            grid on ; 
            V((p-1)*P+m,:) = axis ; 
            xlabel('$j = \log_2 2^j$','Interpreter','Latex','FontSize',fontsize2) ;
            ylabel(['$\log_2 S_{',num2str(p),num2str(p),'}(2^j)$'],'Interpreter','Latex','FontSize',fontsize2)
            pbaspect([1.1,1,1])
            set(gca,'LineWidth',linewidth2)
            xticks(1:4:JMM)
        end
    end
end
clear Vf
Vf(1) = 0 ; Vf(2) = JMM+1 ; Vf(3) = min(V(:,3)) ; Vf(4) = max(V(:,4)) ; 
for p = 1:1:P, subplot(P,P,(p-1)*P+p) ; axis(Vf)  ; end
fig.Position = [57 160 834 643];
sgtitle('Log-wavelet spectrum (diagonal) and wavelet coherence (off-diagonal)','Interpreter','Latex')


%% Univariate analysis
markersize = 15 ; 
linewidth = 2 ; 
fontsize = 20 ;

for j = 1:size(est.WW,3), dWW(j,:) = diag(est.WW(:,:,j)); end
fig = figure(2) ; clf ;   
for p = 1:1:P
    JJ =1:JMM ;
    JJJ = paramsEst.j1:paramsEst.j2 ;
    plot(JJ,log2(abs(dWW(:,p))),'ob-','MarkerSize',markersize,'LineWidth',linewidth); grid on ; hold on ;
    pp = polyfit(JJJ,log2(abs(dWW(JJJ,p)))',1) ; qq = polyval(pp,JJ) ;
    hh(1) = plot(1:JMM,qq,'k--','MarkerSize',markersize,'LineWidth',linewidth) ; 
    hh(2) = plot(JJJ,qq(JJJ),'k-','MarkerSize',markersize,'LineWidth',linewidth) ; 
    xticks(1:2:JMM); xlim([.5 JMM+.5])
end
xlabel('$j = \log_2 2^j$','Interpreter','Latex','FontSize',fontsize) ;
ylabel('$\log_2 S_{m,m}(2^j)$','Interpreter','Latex','FontSize',fontsize)
legend(hh,{'Linear regression','Analysis scales'},'location','best','Interpreter','Latex')
set(gca,'FontSize',fontsize,'LineWidth',linewidth,'TickLabelInterpreter','Latex')

disp(['Univariate H estimates:        [',num2str(diag(est.hU)'-.5),']'])

%% Multivariate analysis
markersize = 15 ; 
linewidth = 2 ; 
fontsize = 20 ;

fig = figure(3) ; clf ; 
set(gca,'FontSize',fontsize) ;   
for p = 1:1:P
    JJ =1:JMM ;
    JJJ = paramsEst.j1:paramsEst.j2 ;
    plot(JJ,log2(estbc.lambdaj(:,p)),'or-','MarkerSize',markersize,'LineWidth',linewidth) ; grid on ; hold on ;
    pp = polyfit(JJJ,log2(estbc.lambdaj(JJJ,p))',1) ;  qq = polyval(pp,JJ) ; 
    hh(1) = plot(1:JMM,qq,'k--','MarkerSize',markersize,'LineWidth',linewidth) ; 
    hh(2) = plot(JJJ,qq(JJJ),'k-','MarkerSize',markersize,'LineWidth',linewidth) ; 
    xticks(1:2:JMM); xlim([.5 JMM+.5])
end
xlabel('$j = \log_2 2^j$','Interpreter','Latex','FontSize',fontsize) ;
ylabel('$\log_2 \lambda_{m}(2^j)$','Interpreter','Latex','FontSize',fontsize)
legend(hh,{'Linear regression','Analysis scales'},'location','best','Interpreter','Latex')
set(gca,'FontSize',fontsize,'LineWidth',linewidth,'TickLabelInterpreter','Latex')

disp(['Multivariate H estimates:   [',num2str(estbc.h),']'])

%% Comparison of analysis
markersize = 12 ; 
linewidth = 2 ; 
fontsize = 20 ;

fig = figure(4) ; clf ;    
for p = 1:1:P
    JJ =1:JMM ;
    JJJ = paramsEst.j1:paramsEst.j2 ;
    hh(1) = plot(JJ,log2(abs(dWW(:,p))),'ob-','MarkerSize',markersize,'LineWidth',linewidth); grid on ; hold on ;
    hh(2) = plot(JJ,log2(estbc.lambdaj(:,p)),'or-','MarkerSize',markersize,'LineWidth',linewidth) ;
    xticks(1:2:JMM); xlim([.5 JMM+.5])
end
xlabel('$j = \log_2 2^j$','Interpreter','Latex','FontSize',fontsize) ;
ylabel('Structure functions','Interpreter','Latex','FontSize',fontsize)
legend(hh,{'$\log_2 S_{m,m}(2^j)$','$\log_2 \lambda_{m}(2^j)$'},'location','best','Interpreter','Latex')
set(gca,'FontSize',fontsize,'LineWidth',linewidth,'TickLabelInterpreter','Latex')