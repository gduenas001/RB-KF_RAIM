
clear; clc; close all;

colors= distinguishable_colors(10);

Nepochs= 100;

load('sigw01v03.mat');
XXerror01= XXerrorsave;
PX01= PXsave;
P_HMI01= P_HMI; 

load('sigw02v03.mat');
XXerror02= XXerrorsave;
PX02= PXsave;
P_HMI02= P_HMI; 

load('sigw025v03.mat');
XXerror025= XXerrorsave;
PX025= PXsave;
P_HMI025= P_HMI; 

load('sigw03v03.mat');
XXerror03= XXerrorsave;
PX03= PXsave;
P_HMI03= P_HMI;

load('sigw04v03.mat');
XXerror04= XXerrorsave;
PX04= PXsave;
P_HMI04= P_HMI;

load('sigwInfv03.mat');
XXerrorInf= XXerrorsave;
PXInf= PXsave;
P_HMIInf= P_HMI; 

load('sigw01v03LS.mat');
XXerrorLS= XXerrorsave;
PXLS= PXsave;
P_HMILS= P_HMI; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the RMSE
RMSE01= sqrt( (1/Nepochs)* sum(XXerror01(1,:).^2) );
RMSE02= sqrt( (1/Nepochs)* sum(XXerror02(1,:).^2) );
RMSE025= sqrt( (1/Nepochs)* sum(XXerror025(1,:).^2) );
RMSE03= sqrt( (1/Nepochs)* sum(XXerror03(1,:).^2) );
RMSE04= sqrt( (1/Nepochs)* sum(XXerror04(1,:).^2) );
RMSEInf= sqrt( (1/Nepochs)* sum(XXerrorInf(1,:).^2) );
RMSELS= sqrt( (1/Nepochs)* sum(XXerrorLS(1,:).^2) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display variances
3*PX01(1,end)
3*PX02(1,end)
3*PX025(1,end)
3*PX03(1,end)
3*PX04(1,end)
3*PXInf(1,end)
3*PXLS(1,end)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on; grid on;
legError01= plot(1:length(XXerrorLS), abs(XXerror01(1,:)), 'b-','linewidth',2);
legError02= plot(1:length(XXerrorLS), abs(XXerror02(1,:)), 'g-','linewidth',2);
legError03= plot(1:length(XXerrorLS), abs(XXerror03(1,:)), 'r-','linewidth',2);
% plot(1:Nepochs, abs(XXerrorInf(1,:)),'k-','linewidth',2);
legErrorLS= plot(1:Nepochs, abs(XXerrorLS(1,:)),'k-','linewidth',2);

% plot(1:0.5:length(PXLS), 3* sqrt(PX01(1,1:end-1)),  '--b','linewidth',2);
% plot(1:0.5:length(PXLS), 3* sqrt(PX02(1,1:end-1)),  '--g','linewidth',2);
% plot(1:0.5:length(PXLS), 3* sqrt(PX03(1,1:end-1)),  '--r','linewidth',2);
% % plot(1:Nepochs, 3* sqrt(PXInf(1,:)),  '--k','linewidth',2);
plot(1:length(PXLS), 3* sqrt(PXLS(1,:)),  '--k','linewidth',2);

% plot(1:Nepochs, -3* sqrt(PX01(1,:)),  '--b','linewidth',2);
% plot(1:Nepochs, -3* sqrt(PX02(1,:)),  '--g','linewidth',2);
% plot(1:Nepochs, -3* sqrt(PX03(1,:)),  '--r','linewidth',2);
% % plot(1:Nepochs, -3* sqrt(PXInf(1,:)),  '--k','linewidth',2);
% plot(1:Nepochs, -3* sqrt(PXLS(1,:)),  '--k','linewidth',2);

set(gca,'fontsize',8)

xlabel({'time epochs'},'interpreter','latex','fontsize',10);
ylabel({'$\hat{\epsilon}$'},'interpreter','latex','fontsize',10);
legend([legError01,legError02,legError03,legErrorLS],...
    {'$\sigma_w = 0.1$','$\sigma_w = 0.2$','$\sigma_w = 0.3$','LS or $\sigma_w = \infty$'},...
    'location','northeast', 'interpreter','latex','fontsize',10);

% % Save figure
% fig = gcf;
% fig.PaperUnits = 'centimeters';
% fig.PaperPosition = [0 0 9 7];
% print('error&variance','-dpdf','-r0')
%}




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on; grid on;
set(gca,'fontsize',8);
xlabel({'time epochs'},'interpreter','latex','fontsize',10); 
ylabel({'P(HMI)'},'interpreter','latex','fontsize',10); 
set(gca,'Yscale','log');
ylim([1e-10,1e-4]);
legHMI01= plot(1:Nepochs, P_HMI01,'-', 'color', colors(1,:), 'linewidth',2);
legHMI02= plot(1:Nepochs, P_HMI02,'-', 'color', colors(2,:), 'linewidth',2);
% legHMI025= plot(1:Nepochs, P_HMI025,'-', 'color', colors(3,:), 'linewidth',2);
legHMI03= plot(1:Nepochs, P_HMI03,'-', 'color', colors(4,:), 'linewidth',2);
legHMI04= plot(1:Nepochs, P_HMI04,'-', 'color', colors(5,:), 'linewidth',2);
legHMIInf= plot(1:Nepochs, P_HMIInf,'-', 'color', colors(6,:), 'linewidth',2);
legHMILS= plot(1:Nepochs, P_HMILS,'--', 'color', colors(7,:), 'linewidth',2);


% legend([legHMI01,legHMI02,legHMI03,legHMI04,legHMIInf,legHMILS],...
%     {'$\sigma_w = 0.1$','$\sigma_w = 0.2$','$\sigma_w = 0.3$','$\sigma_w = 0.4$',...
%     '$\sigma_w = \infty$','Snapshot'},...
%     'location','northeast', 'interpreter','latex','fontsize',10);

% % Save figure
% fig = gcf;
% fig.PaperUnits = 'centimeters';
% fig.PaperPosition = [0 0 9 7];
% print('P_HMI','-dpdf','-r0')
% %}

