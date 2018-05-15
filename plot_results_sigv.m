
clear; clc; close all;

Nepochs= 100;

load('sigw03v01.mat');
XXerror01= XXerrorsave;
PX01= PXsave;
P_HMI01= P_HMI; 

load('sigw03v02.mat');
XXerror02= XXerrorsave;
PX02= PXsave;
P_HMI02= P_HMI; 

load('sigw03v025.mat');
XXerror025= XXerrorsave;
PX025= PXsave;
P_HMI025= P_HMI; 

load('sigw03v03.mat');
XXerror03= XXerrorsave;
PX03= PXsave;
P_HMI03= P_HMI; 

load('sigw03v04.mat');
XXerror04= XXerrorsave;
PX04= PXsave;
P_HMI04= P_HMI; 

load('sigw01v01LS.mat');
XXerror01LS= XXerrorsave;
PX01LS= PXsave;
P_HMI01LS= P_HMI; 

load('sigw01v02LS.mat');
XXerror02LS= XXerrorsave;
PX02LS= PXsave;
P_HMI02LS= P_HMI; 

load('sigw01v025LS.mat');
XXerror025LS= XXerrorsave;
PX025LS= PXsave;
P_HMI025LS= P_HMI; 

load('sigw01v03LS.mat');
XXerror03LS= XXerrorsave;
PX03LS= PXsave;
P_HMI03LS= P_HMI; 

load('sigw01v04LS.mat');
XXerror04LS= XXerrorsave;
PX04LS= PXsave;
P_HMI04LS= P_HMI; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the RMSE
RMSE01= sqrt( (1/Nepochs)* sum(XXerror01(1,:).^2) );
RMSE02= sqrt( (1/Nepochs)* sum(XXerror02(1,:).^2) );
RMSE025= sqrt( (1/Nepochs)* sum(XXerror025(1,:).^2) );
RMSE03= sqrt( (1/Nepochs)* sum(XXerror03(1,:).^2) );
RMSE04= sqrt( (1/Nepochs)* sum(XXerror04(1,:).^2) );

RMSE01LS= sqrt( (1/Nepochs)* sum(XXerror01LS(1,:).^2) );
RMSE02LS= sqrt( (1/Nepochs)* sum(XXerror02LS(1,:).^2) );
RMSE025LS= sqrt( (1/Nepochs)* sum(XXerror025LS(1,:).^2) );
RMSE03LS= sqrt( (1/Nepochs)* sum(XXerror03LS(1,:).^2) );
RMSE04LS= sqrt( (1/Nepochs)* sum(XXerror04LS(1,:).^2) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display the variances
3*PX01(1,end)
3*PX02(1,end)
3*PX025(1,end)
3*PX03(1,end)
3*PX04(1,end)

3*PX01LS(1,end)
3*PX02LS(1,end)
3*PX025LS(1,end)
3*PX03LS(1,end)
3*PX04LS(1,end)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
figure; hold on; grid on;
legError01= plot(1:Nepochs, abs(XXerror01(1,:)), 'b-','linewidth',2);
legError02= plot(1:Nepochs, abs(XXerror02(1,:)), 'r-','linewidth',2);
legError025= plot(1:Nepochs, abs(XXerror025(1,:)),'k-','linewidth',2);
legError03= plot(1:Nepochs, abs(XXerror03(1,:)), 'm-','linewidth',2);
legError04= plot(1:Nepochs, abs(XXerror04(1,:)), 'g-','linewidth',2);
legError01LS= plot(1:Nepochs, abs(XXerror01LS(1,:)),'b-*','linewidth',2);
legError02LS= plot(1:Nepochs, abs(XXerror02LS(1,:)),'r-*','linewidth',2);
legError025LS= plot(1:Nepochs, abs(XXerror025LS(1,:)),'k-*','linewidth',2);
legError03LS= plot(1:Nepochs, abs(XXerror03LS(1,:)),'m-*','linewidth',2);
legError04LS= plot(1:Nepochs, abs(XXerror04LS(1,:)),'g-*','linewidth',2);

plot(1:Nepochs, 3* sqrt(PX01(1,:)),  '--b','linewidth',2);
plot(1:Nepochs, 3* sqrt(PX02(1,:)),  '--r','linewidth',2);
plot(1:Nepochs, 3* sqrt(PX025(1,:)),  '--k','linewidth',2);
plot(1:Nepochs, 3* sqrt(PX03(1,:)),  '--m','linewidth',2);
plot(1:Nepochs, 3* sqrt(PX04(1,:)),  '--g','linewidth',2);
plot(1:Nepochs, 3* sqrt(PX01LS(1,:)),  '--*b','linewidth',2);
plot(1:Nepochs, 3* sqrt(PX02LS(1,:)),  '--*r','linewidth',2);
plot(1:Nepochs, 3* sqrt(PX025LS(1,:)),  '--*k','linewidth',2);
plot(1:Nepochs, 3* sqrt(PX03LS(1,:)),  '--*m','linewidth',2);
plot(1:Nepochs, 3* sqrt(PX04LS(1,:)),  '--*g','linewidth',2);
set(gca,'fontsize',8)

xlabel({'time epochs'},'interpreter','latex','fontsize',10);
ylabel({'$\hat{\epsilon}$'},'interpreter','latex','fontsize',10);
legend([legError01,legError02,legError025,legError03,legError04,...
    legError01LS,legError02LS,legError025LS,legError03LS,legError04LS],...
    {'KF with $\sigma_v = 0.1$','KF with $\sigma_v = 0.2$','KF with $\sigma_v = 0.25$','KF with $\sigma_v = 0.3$','KF with $\sigma_v = 0.4$',...
    'LS with $\sigma_v = 01$','LS with $\sigma_v = 02$','LS with $\sigma_v = 0.25$','LS with $\sigma_v = 0.3$','LS with $\sigma_v = 0.4$'},....
    'location','northeast', 'interpreter','latex','fontsize',10);


% Save figure
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 9 7];
print('error&variance','-dpdf','-r0')
%}


figure; hold on; grid on;
set(gca,'fontsize',8);
xlabel({'time epochs'},'interpreter','latex','fontsize',10); 
ylabel({'P(HMI)'},'interpreter','latex','fontsize',10); 
set(gca,'Yscale','log');
% ylim([1e-10,1e-4]);
% legHMI01= plot(1:Nepochs, P_HMI01,'-b', 'linewidth',2);
legHMI02= plot(1:Nepochs, P_HMI02,'-b', 'linewidth',2);
legHMI025= plot(1:Nepochs, P_HMI025,'-k', 'linewidth',2);
legHMI03= plot(1:Nepochs, P_HMI03,'-r', 'linewidth',2);
legHMI04= plot(1:Nepochs, P_HMI04,'-g', 'linewidth',2);
legHMI01LS= plot(1:Nepochs, P_HMI01LS,'--b', 'linewidth',2);
% legHMI02LS= plot(1:Nepochs, P_HMI02LS,'--r', 'linewidth',2);
legHMI025LS= plot(1:Nepochs, P_HMI025LS,'--k', 'linewidth',2);
legHMI03LS= plot(1:Nepochs, P_HMI03LS,'--r', 'linewidth',2);
legHMI04LS= plot(1:Nepochs, P_HMI04LS,'--g', 'linewidth',2);

% legend([legHMI01,legHMI025,legHMI03,legHMI04,...
%         legHMI01LS,legHMI025LS,legHMI03LS,legHMI04LS],...
%     {'$\sigma_v = 0.1, 0.2$','$\sigma_v = 0.25$','$\sigma_v = 0.3$','$\sigma_v = 0.4$',...
%     'LS with $\sigma_v = 0.1, 0.2$','LS with $\sigma_v = 0.25$','LS with $\sigma_v = 0.3$','LS with $\sigma_v = 0.4$',},...
%     'location','northeast', 'interpreter','latex','fontsize',10);
legend([legHMI02,legHMI025,legHMI03,legHMI04],...
    {'$\sigma_v = 0.2$','$\sigma_v = 0.25$','$\sigma_v = 0.3$','$\sigma_v = 0.4$'},...
    'location','northeast', 'interpreter','latex','fontsize',10);

% % Save figure
% fig = gcf;
% fig.PaperUnits = 'centimeters';
% % fig.PaperPosition = [0 0 9 7];
% fig.PaperPosition = [0 0 9 3];
% print('P_HMI_sig_v_zoom','-dpdf','-r0')
% %}

