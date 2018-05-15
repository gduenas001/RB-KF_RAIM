
clear; clc; close all;

% input parameters
sig_w_real= 0.1;
sig_w= 0.2;
sig_v= 0.3;
alert_limit= 1;
Pi_fail= 5e-4;
vel_max= 40; vel_max= vel_max/3.6; % km/h to m/s
dt= 0.1;
n_F= 6;

% input paramters that do not affect performance
Nepochs= 100;
m= 2;
m_F= 2;
C_REQ= 1e-7;
alpha= [1; 0];
I_H= 1e-10;

% build paramters
n= n_F * m_F;
W_real= kron(eye(m),sig_w_real^2);
W= kron(eye(m),sig_w^2);
V= kron(eye(n),sig_v^2);
T_RB= chi2inv(1 - C_REQ,n);
rng(1); % random seed for repeatability

% Generate landmarks
LM= [(rand(1,30)-0.5) * 20;
      rand(1,30) * 40];
LM= LM(:,1:n_F);

% n_max
for r=1:n_F
    if nchoosek(n_F,r) * Pi_fail^r < I_H
        n_max= r-1;
        break;
    end
    if r > n_F*m_F - m
        error('We cannot monitor all modes of failure');
    end
end
  
% Initialization
xtrue= [0;0];
XX= xtrue;
PX= kron(eye(m),eps);


PH0p= 1; % No previous faults probability
xtruesave= ones(m,Nepochs);
XXerrorsave= ones(m,Nepochs) * (-1);
XXsave= ones(m,Nepochs); XXsave(:,1)= XX;
PXsave= ones(m,Nepochs); PXsave(:,1)= diag(PX);
P_HMI= ones(1,Nepochs);
for epoch= 1:Nepochs
    %% prediction step
    % Increase speed in y at every epoch
    vel= [vel_max * cosd( epoch * 360/100 ) ;
          vel_max * sind( epoch * 360/100 )];
    
    % Generate noise and predict
    w= mvnrnd(zeros(m,1), W_real)';
    xtrue= xtrue + vel*dt + w;
    xtruesave(:,epoch)= xtrue;
    XX= XX + vel*dt;
    PX= PX + W;
        
    %% update step
    % Build matrices
    H= repmat([-1,0;0,-1],n_F,1);
    D= [H; eye(m)];
    Delta= [V,          zeros(n,m);
            zeros(m,n), PX];
    invDelta= inv(Delta);
    S= (D'*invDelta*D) \ D' * invDelta;
    DeltaI_DS= invDelta - invDelta * D * S;
    
    % obtain msmts
    v= mvnrnd(zeros(n,1), V)';
    z= H * xtrue + LM(:) + v;

    % Composed measurement
    z_star= z - LM(:);
    l= [z_star; XX];
    
    %% RAIM
    
    % Important parameters
    PX_hat= S * Delta * S';
    sig_hat= sqrt(alpha'*PX_hat*alpha);
    P_HMI_H_fn_generic= @(f_norm,fhat_eps,lambda2) (-1)*...
        ( 1 - cdf('normal', alert_limit, fhat_eps*f_norm, sig_hat) + ...
        cdf('normal', -alert_limit, fhat_eps*f_norm, sig_hat) ) * ...
        cdf('Noncentral Chi-square', T_RB, n, lambda2*f_norm^2);
    
    % Fault probabilities
    PH0= (1 - Pi_fail)^n_F;
    PH= ones(1,n_max) * (-1);
    P_HMI_H0p= ones(1,n_max) * (-1);
    P_HMI_H1p= ones(1,n_max) * (-1);
    for n_f= 1:n_max
        % Probability of current fault with n_f faulted msmts
        PH(n_f)= nchoosek(n_F,n_f) * PH0 * ( Pi_fail / (1-Pi_fail) )^n_f;
        
        %% H0p ------ No previous faults
        f_l_u= worst_case_fault(n_F,m_F,m,n_f,false, DeltaI_DS,S,alpha);
        f_l_u= f_l_u / norm(f_l_u); % normalize -- I don't know if necessary
        
        % Worst case fault norms
        fhat_eps= - alpha' * S * f_l_u;
        lambda2= f_l_u' * DeltaI_DS * f_l_u;
        P_HMI_H_fn= @(f_norm) P_HMI_H_fn_generic(f_norm,fhat_eps,lambda2);
        [f_norm, P_HMI_H]= fminbnd(P_HMI_H_fn, 0, 10);
        P_HMI_H0p(n_f)= - P_HMI_H;
        
        %% H1p ------ Yes previous faults
        f_l_u= worst_case_fault(n_F,m_F,m,n_f,true, DeltaI_DS,S,alpha);
        f_l_u= f_l_u / norm(f_l_u); % normalize -- I don't know if necessary
        
        % Worst case fault norms
        fhat_eps= - alpha' * S * f_l_u;
        lambda2= f_l_u' * DeltaI_DS * f_l_u;
        P_HMI_H_fn= @(f_norm) P_HMI_H_fn_generic(f_norm,fhat_eps,lambda2);
        [f_norm, P_HMI_H]= fminbnd(P_HMI_H_fn, 0, 10);
        P_HMI_H1p(n_f)= - P_HMI_H;
    end
    
    % P(H0,H1p) ------------------- not current faults
    E_H0H1p= zeros(m,n+m);
    E_H0H1p(end-(m-1):end,end-(m-1):end)= eye(m);
    f_l_u= E_H0H1p' * ( (E_H0H1p * (DeltaI_DS) * E_H0H1p' ) \ E_H0H1p * S' * alpha);
    f_l_u= f_l_u / norm(f_l_u);
    
    % Worst case fault norms
    fhat_eps= - alpha' * S * f_l_u;
    lambda2= f_l_u' * DeltaI_DS * f_l_u;
    P_HMI_H_fn= @(f_norm) P_HMI_H_fn_generic(f_norm,fhat_eps,lambda2);
    [f_norm, P_HMI_H]= fminbnd(P_HMI_H_fn, 0, 10);
    P_HMI_H0H1p= - P_HMI_H;
    
    % P(H0,H0p) ------------------- not faults ever
    P_HMI_H0H0p= 2 * cdf('normal', -alert_limit, 0, sig_hat) *...
        cdf('Noncentral Chi-square', T_RB, n, 0);
    
    % Total P(HMI) ------------------- 
    P_HMI(epoch)= PH0p * (PH0*P_HMI_H0H0p + sum(PH .* P_HMI_H0p)) + ...
                  (1 - PH0p) * (PH0*P_HMI_H0H1p + sum(PH .* P_HMI_H1p)) + ...
                  I_H; % nchoosek(n_F,n_max+1) * Pi_fail^(n_max+1);   
    %% Update
    XX= S * l;
    PX= S * Delta * S';
    
    % Update prob of no failures
    PH0p= PH0p * (1 - Pi_fail)^n;

    % Save values
    XXerrorsave(:,epoch)= XX - xtrue;
    XXsave(:,epoch)= XX;
    PXsave(:,epoch)= diag(PX);
    
end




%% Plots



% X - Y positions and landmarks
% subplot(2,3,3); hold on; grid on;
figure; hold on;
xlabel('X [m]'); ylabel('Y [m]');
plot(xtruesave(1,:),xtruesave(2,:),'g-');
% plot(XXsave(1,:),XXsave(2,:),'b-');
plot(LM(1,:),LM(2,:),'k+','markersize',10);
axis equal
% legend({'$\textbf{x}$','$\hat{\textbf{x}}$','landmarks'},...
%     'location','northeast', 'interpreter','latex');
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 6 6];
print('lm_map','-dpdf','-r0')
%}

%{
% X-Y variance envelopes
subplot(2,3,[4,5,6]); hold on; grid on;
xlabel('time epochs'); ylabel('m');
% X - coordinate
leg.x= plot(1:Nepochs, xtruesave(1,:),'-g','linewidth',2);
leg.xhat= plot(1:Nepochs, XXsave(1,:),'-b');
leg.sigx= plot(1:Nepochs, XXsave(1,:) + 3*sqrt(PXsave(1,:)),'--b');
plot(1:Nepochs, XXsave(1,:) - 3*sqrt(PXsave(1,:)),'--b');
% Y - coordinate
leg.y= plot(1:Nepochs, xtruesave(2,:),'-r','linewidth',2);
leg.yhat= plot(1:Nepochs, XXsave(2,:),'-b');
plot(1:Nepochs, XXsave(2,:) + 3*sqrt(PXsave(2,:)),'--b');
plot(1:Nepochs, XXsave(2,:) - 3*sqrt(PXsave(2,:)),'--b');
legend([leg.x,leg.xhat,leg.sigx,leg.y,leg.yhat],...
    {'$x$','$\hat{x}$','$3 \hat{\sigma}$','$y$','$\hat{y}$'},...
    'location','northwest', 'interpreter','latex');
%}

% X-Y standard deviations
% subplot(2,3,[4,5,6]); hold on; grid on;
fig_SD= figure; hold on; grid on;
xlabel('time epochs'); ylabel('m');
legKF.errorx= plot(1:Nepochs, abs(XXerrorsave(1,:)),'b-','linewidth',2);
legKF.sigx= plot(1:Nepochs, 3* sqrt(PXsave(1,:)),'--b','linewidth',2);
legend([legKF.errorx,legKF.sigx],...
    {'$\hat{\epsilon}_x$','$3 \hat{\sigma}_x$'},...
    'location','northwest', 'interpreter','latex');


% P(HMI)
% subplot(2,3,[1,2]); hold on; grid on;
fig_HMI= figure; hold on; grid on;
xlabel('time epochs'); ylabel('P(HMI)');
set(gca,'Yscale','log');
plot(1:Nepochs, P_HMI,'-b','linewidth',2);
ylim([1e-10,1e-4]);
% title(['$n_F$ = ',num2str(n_F), '  $\sigma_w$ = ',num2str(sig_w), '  $\sigma_v$ = ',num2str(sig_v)],...
%     'interpreter','latex');

% Make it fullscreen
% set(gcf, 'Position', get(0, 'Screensize'));


% % Save figure
% name= 'KF_vs_LS';
% saveas(gcf,name);
% saveas(gcf,[name,'.png']);
% saveas(gcf,['P_HMI_sig_w-',num2str()]);
% saveas(gcf,['P_HMI_sig_w-',num2str(sig_w),'.png']);


% save('sigw02v03','P_HMI','PXsave','XXerrorsave');

