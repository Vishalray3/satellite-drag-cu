%% plots BFF
clc
clearvars
load Truth_cube_dria_iner_bff
load BFF_orbit_Kfirst
Cd1 = Cd_true;
Af1 = Af_analytical;
load Truth_cube_dria_iner_bff_Klast
load BFF_orbit_Klast
Cd2 = Cd_true;
Af2 = Af_analytical;
%%
n_plot = ceil(T_orb/del_T);
xvec = 0:360/n_plot:360-360/n_plot;
figure(1)
plot(xvec, Cd1(1:n_plot), 'k','LineWidth', 1)
hold on
plot(xvec, Cd2(1:n_plot), '--r','LineWidth', 1)
xlabel('Eccentric anomaly (deg)')
ylabel('Drag-coefficient')
legend('K = 1e6', 'K = 1e8')
title('Drag-coefficient in orbit')
set(gca,'FontSize',18)

%%
figure(2)
subplot(3,1,1)
plot( Af1(1,:), 'k','LineWidth', 1)
hold on
plot( Af2(1,:), '--r','LineWidth', 1)
ylabel('$\overline{\mathcal{A}}_0$','interpreter','latex')
title('Fourier coefficients')
legend('K = 1e6', 'K = 1e8')
set(gca,'FontSize',18)

subplot(3,1,2)
plot( Af1(5,:), 'k','LineWidth', 1)
hold on
plot( Af2(5,:), '--r','LineWidth', 1)
ylabel('$\overline{\mathcal{A}}_4$','interpreter','latex')
set(gca,'FontSize',18)

subplot(3,1,3)
plot( Af1(9,:), 'k','LineWidth', 1)
hold on
plot( Af2(9,:), '--r','LineWidth', 1)
ylabel('$\overline{\mathcal{A}}_8$','interpreter','latex')
xlabel('Eccentric anomaly (deg)')
set(gca,'FontSize',18)
%%
xvec = time_prop_utc/T_orb;
% load Case3_gpm_full_iner
figure(1)
plot(xvec, obs_rank(1,:), 'k','LineWidth', 1)
xlabel('Orbits')
ylabel('Rank')
title('Rank of observability matrix')
grid on
set(gca,'FontSize',18)
% set(gca,'XTick',xvec(1):2:xvec(end))
% ylim([1 15])

figure(2)
semilogy(xvec, s_obs(end,:), 'k','LineWidth', 1)
hold on
semilogy(xvec, s_obs(end-1,:), '--r','LineWidth', 1)
semilogy(xvec, s_thresh, '-.g','LineWidth', 1)
leg = legend('$s_1$','$s_2$','threshold');
set(leg,'interpreter','latex')
xlabel('Orbits')
grid on
ylabel('Singular values')
title('Singular values of observability matrix')
set(gca,'FontSize',18)
set(gca,'XTick',xvec(1):2:xvec(end))
%% standard deviation
% xvec = time_prop_utc/T_orb;
% figure(3)
% ax1 = subplot(2,1,1);
% plot(xvec, sigma_x(7,:)/abs(X_nom(7)), 'k','LineWidth', 1)
% hold on
% semilogy(xvec, sigma_x(9,:)/abs(X_nom(9)), 'r','LineWidth', 1)
% semilogy(xvec, sigma_x(11,:)/abs(X_nom(11)), 'b','LineWidth', 1)
% semilogy(xvec, sigma_x(29,:)/abs(X_nom(29)), 'g','LineWidth', 1)
% title('Uncertainty of Fourier coefficients')
% leg = legend('$\overline{\mathcal{A}}_0$','$\overline{\mathcal{A}}_2$','$\overline{\mathcal{A}}_4$','$\overline{\mathcal{B}}_1$');
% set(leg,'interpreter','latex');
% set(gca,'FontSize',18)
% % set(gca,'YTick',[1e-3,1e-2,1e-1])
% grid on
% ax2 = subplot(2,1,2);
% plot(xvec, sigma_x(8,:)/abs(X_nom(8)), 'k','LineWidth', 1)
% hold on
% plot(xvec, sigma_x(10,:)/abs(X_nom(10)), 'r','LineWidth', 1)
% plot(xvec, sigma_x(13,:)/abs(X_nom(13)), 'b','LineWidth', 1)
% leg = legend('$\overline{\mathcal{A}}_1$','$\overline{\mathcal{A}}_{3}$','$\overline{\mathcal{A}}_{6}$');
% set(leg,'interpreter','latex');
% xlabel('Orbits')
% set(gca,'FontSize',18)
% % ylim([0 0.2])
% grid on
% p1=get(ax1,'position');
% p2=get(ax2,'position');
% height=p1(2)+p1(4)-p2(2);
% h3=axes('position',[p2(1) p2(2) p2(3) height],'visible','off');
% h_label=ylabel('Standard deviation','visible','on');
% set(gca,'FontSize',18)

xvec = time_prop_utc/T_orb;
figure(3)
ax1 = subplot(2,1,1);
plot(xvec, sigma_x(7,:)/abs(X_nom(7)), 'k','LineWidth', 1)
hold on
semilogy(xvec, sigma_x(8,:)/abs(X_nom(8)), 'r','LineWidth', 1)
title('Uncertainty of Fourier coefficients')
leg = legend('$\overline{\mathcal{A}}_0$','$\overline{\mathcal{A}}_4$');
set(leg,'interpreter','latex');
set(gca,'FontSize',18)
% set(gca,'YTick',[1e-3,1e-2,1e-1])
grid on
ax2 = subplot(2,1,2);
plot(xvec, sigma_x(9,:)/abs(X_nom(9)), 'k','LineWidth', 1)
hold on
plot(xvec, sigma_x(10,:)/abs(X_nom(10)), 'r','LineWidth', 1)
plot(xvec, sigma_x(11,:)/abs(X_nom(11)), 'g','LineWidth', 1)
leg = legend('$\overline{\mathcal{A}}_8$','$\overline{\mathcal{A}}_{12}$','$\overline{\mathcal{A}}_{16}$');
set(leg,'interpreter','latex');
xlabel('Orbits')
set(gca,'FontSize',18)
% ylim([0 0.2])
grid on
p1=get(ax1,'position');
p2=get(ax2,'position');
height=p1(2)+p1(4)-p2(2);
h3=axes('position',[p2(1) p2(2) p2(3) height],'visible','off');
h_label=ylabel('Standard deviation','visible','on');
set(gca,'FontSize',18)


%% Correlation matrix
figure()
imagesc(abs(Corr_mat([7,9,11,13,28,30,68,92,97],[7,9,11,13,28,30,68,92,97])))
colormap copper
c = colorbar;
c.Label.String = 'Correlation coefficients';
title('Correlation matrix')
xticks([1:9])
yticks([1:9])
set(gca,'TickLabelInterpreter', 'latex');
% cube
% set(gca,'XTickLabel', {'$\overline{\mathcal{A}}_0$','$\overline{\mathcal{A}}_4$','$\overline{\mathcal{A}}_8$','$\overline{\mathcal{A}}_{12}$',...
%     '$\overline{\mathcal{A}}_{16}$','$\overline{\mathcal{A}}_{20}$','$\overline{\mathcal{A}}_{24}$','$\overline{\mathcal{A}}_{28}$'});
% set(gca,'YTickLabel', {'$\overline{\mathcal{A}}_0$','$\overline{\mathcal{A}}_4$','$\overline{\mathcal{A}}_8$','$\overline{\mathcal{A}}_{12}$',...
%     '$\overline{\mathcal{A}}_{16}$','$\overline{\mathcal{A}}_{20}$','$\overline{\mathcal{A}}_{24}$','$\overline{\mathcal{A}}_{28}$'});

%gpm 
% set(gca,'XTickLabel', {'$\overline{\mathcal{A}}_0$','$\overline{\mathcal{A}}_1$','$\overline{\mathcal{A}}_2$','$\overline{\mathcal{A}}_{3}$',...
%     '$\overline{\mathcal{A}}_{4}$','$\overline{\mathcal{A}}_{6}$','$\overline{\mathcal{B}}_{1}$'});
% set(gca,'YTickLabel', {'$\overline{\mathcal{A}}_0$','$\overline{\mathcal{A}}_1$','$\overline{\mathcal{A}}_2$','$\overline{\mathcal{A}}_{3}$',...
%     '$\overline{\mathcal{A}}_{4}$','$\overline{\mathcal{A}}_{6}$','$\overline{\mathcal{B}}_{1}$'});
% gpm bodf
set(gca,'XTickLabel', {'$\overline{A}_{0,0}$','$\overline{A}_{2,0}$','$\overline{A}_{4,0}$','$\overline{A}_{6,0}$','$\overline{A}_{0,1}$','$\overline{A}_{2,1}$',...
   '$\overline{A}_{2,3}$','$\overline{B}_{1,0}$','$\overline{B}_{1,1}$'});
set(gca,'YTickLabel', {'$\overline{A}_{0,0}$','$\overline{A}_{2,0}$','$\overline{A}_{4,0}$','$\overline{A}_{6,0}$','$\overline{A}_{0,1}$','$\overline{A}_{2,1}$',...
   '$\overline{A}_{2,3}$','$\overline{B}_{1,0}$','$\overline{B}_{1,1}$'});
set(gca,'FontSize',18)

%% consider covariance results

load bff_iner_consider_order0
sig_pos_est = [norm(sigma_x(1:3,end)); norm(sigma_xx(1:3,end))];
sig_vel_est = [norm(sigma_x(4:6,end)); norm(sigma_xx(4:6,end))];
sig_cd_est = [norm(sigma_x(7,end)); norm(sigma_xx(7,end))];
load bff_iner_consider_order4
sig_pos_est = [sig_pos_est,[norm(sigma_x(1:3,end)); norm(sigma_xx(1:3,end))] ];
sig_vel_est = [sig_vel_est,[norm(sigma_x(4:6,end)); norm(sigma_xx(4:6,end))]];
sig_cd_est = [sig_cd_est,[norm(sigma_x(7,end));norm(sigma_xx(7,end))]];
load bff_iner_consider_order8
sig_pos_est = [sig_pos_est,[norm(sigma_x(1:3,end)); norm(sigma_xx(1:3,end))] ];
sig_vel_est = [sig_vel_est,[norm(sigma_x(4:6,end)); norm(sigma_xx(4:6,end))]];
sig_cd_est = [sig_cd_est,[norm(sigma_x(7,end));norm(sigma_xx(7,end))]];
load bff_iner_consider_order12
sig_pos_est = [sig_pos_est,[norm(sigma_x(1:3,end)); norm(sigma_xx(1:3,end))] ];
sig_vel_est = [sig_vel_est,[norm(sigma_x(4:6,end)); norm(sigma_xx(4:6,end))]];
sig_cd_est = [sig_cd_est,[norm(sigma_x(7,end));norm(sigma_xx(7,end))]];
% load Case1_sphere_consider_order4_cbar
% sig_pos_est = [sig_pos_est,[norm(sigma_x(1:3,end)); norm(sigma_xx(1:3,end))] ];
% sig_vel_est = [sig_vel_est,[norm(sigma_x(4:6,end)); norm(sigma_xx(4:6,end))]];
% sig_cd_est = [sig_cd_est,[norm(sigma_x(7,end));norm(sigma_xx(7,end))]];

vec = [0,4,8,12];
figure()
semilogy(vec, sig_pos_est(1,:),'-k','LineWidth',1)
hold on
semilogy(vec, sig_vel_est(1,:),'-r','LineWidth',1)
semilogy(vec, sig_cd_est(1,:),'-g','LineWidth',1)
semilogy(vec, sig_pos_est(2,:),'--k','LineWidth',1)
semilogy(vec, sig_vel_est(2,:),'--r','LineWidth',1)
semilogy(vec, sig_cd_est(2,:),'--g','LineWidth',1)
legend({'Position','Velocity','Order 0'})
xlabel('Truncation order for estimation')
ylabel('Standard deviation')
title('Estimated and considered uncertainties')
set(gca,'FontSize',18)
set(gca,'XTick',vec)


 %% Monte Carlo plots 
% cases to denote for the paper: {A0} (Case 1), {A0, A4} (Case 2),
% {A0,A4,A8} (Case 3), {A0,A4,A8,A12} (Case 5), {A0,A8} (Case 4)
str_case = 'MC_bff_case';
ind_case = [1,2,3,5,4];
for ii = 1:numel(ind_case)
    load(strcat(str_case,num2str(ii),'_mod'))
    cd_mean(1,ind_case(ii)) = mean(mc_cderr);
    cd_std(1,ind_case(ii)) = std(mc_cderr);
    pos_mean(1,ind_case(ii)) = mean(mc_poserr);
    pos_std(1,ind_case(ii)) = std(mc_poserr);   
    vel_mean(1,ind_case(ii)) = mean(mc_velerr);
    vel_std(1,ind_case(ii)) = std(mc_velerr);   
    rms_ratio_mean(1,ind_case(ii)) = mean(mc_rms_res_norm);
    rms_ratio_std(1,ind_case(ii)) = std(mc_rms_res_norm);
    
    load(strcat(str_case,num2str(ii)))
    cd_mean(2,ind_case(ii)) = mean(mc_cderr);
    cd_std(2,ind_case(ii)) = std(mc_cderr);
    pos_mean(2,ind_case(ii)) = mean(mc_poserr);
    pos_std(2,ind_case(ii)) = std(mc_poserr);   
    vel_mean(2,ind_case(ii)) = mean(mc_velerr);
    vel_std(2,ind_case(ii)) = std(mc_velerr);   
    rms_ratio_mean(2,ind_case(ii)) = mean(mc_rms_res_norm);
    rms_ratio_std(2,ind_case(ii)) = std(mc_rms_res_norm); 
    
    load(strcat(str_case,num2str(ii), '_mod_noErrOFF'))
    cd_mean(3,ind_case(ii)) = mean(mc_cderr);
    cd_std(3,ind_case(ii)) = std(mc_cderr);
    pos_mean(3,ind_case(ii)) = mean(mc_poserr);
    pos_std(3,ind_case(ii)) = std(mc_poserr);   
    vel_mean(3,ind_case(ii)) = mean(mc_velerr);
    vel_std(3,ind_case(ii)) = std(mc_velerr);   
    rms_ratio_mean(3,ind_case(ii)) = mean(mc_rms_res_norm);
    rms_ratio_std(3,ind_case(ii)) = std(mc_rms_res_norm);      
end
 
figure(1)
errorbar(cd_mean(1,:),cd_std(1,:),'k*','LineWidth',3)
hold on
errorbar(cd_mean(2,:),cd_std(2,:),'ro','LineWidth',3)
errorbar(cd_mean(3,:),cd_std(3,:),'gx','LineWidth',3)
ylabel('Drag-coefficient error')
xlabel('Case')
% set(gca,'XTickLabel', {'Case 1','Case 2','Case 3','Case 4','Case 5','Case 6','Case 7'});
title('Drag-coefficient error')
legend('Modeled higher-orders','Ignored higher-orders','True OFF coefficients')
set(gca,'FontSize',18)
grid on

% axes('position',[.3 .4 .47 .25])
% box on 
% errorbar(2:6,cd_mean(1,2:6),cd_std(1,2:6),'k*','LineWidth',1)
% % hold on
% % errorbar(2:6,cd_mean(2,2:6),cd_std(2,2:6),'ro','LineWidth',1)
% set(gca,'FontSize',16)
% % axis tight
% grid on

figure(2)
errorbar(pos_mean(1,:),pos_std(1,:),'k*','LineWidth',3)
hold on
errorbar(pos_mean(2,:),pos_std(2,:),'ro','LineWidth',3)
errorbar(pos_mean(3,:),pos_std(3,:),'gx','LineWidth',3)
ylabel('Initial position error (m)')
% set(gca,'XTickLabel', {'Case 1','Case 2','Case 3','Case 4','Case 5','Case 6','Case 7'});
title('Initial position error')
legend('Modeled higher-orders','Ignored higher-orders','True OFF coefficients')
set(gca,'FontSize',18)
grid on

% axes('position',[.3 .4 .5 .25])
% box on 
% errorbar(2:7,pos_mean(1,2:7),pos_std(1,2:7),'k*','LineWidth',1)
% % hold on
% % errorbar(2:7,pos_mean(2,2:7),pos_std(2,2:7),'ro','LineWidth',1)
% set(gca,'FontSize',16)
% % axis tight
% grid on
% xticks([2:7])

figure(3)
errorbar(rms_ratio_mean(1,:),rms_ratio_std(1,:),'k*','LineWidth',3)
hold on
errorbar(rms_ratio_mean(2,:),rms_ratio_std(2,:),'ro','LineWidth',3)
errorbar(rms_ratio_mean(3,:),rms_ratio_std(3,:),'gx','LineWidth',3)
ylabel('Measurement residual ratio')
% set(gca,'XTickLabel', {'Case 1','Case 2','Case 3','Case 4','Case 5','Case 6','Case 7'});
title('Measurement residual ratio')
legend('Modeled higher-orders','Ignored higher-orders','True OFF coefficients')
set(gca,'FontSize',18)
grid on

axes('position',[.3 .4 .5 .25])
box on 
errorbar(2:4,rms_ratio_mean(1,2:4),rms_ratio_std(1,2:4),'k*','LineWidth',3)
hold on
errorbar(2:4,rms_ratio_mean(2,2:4),rms_ratio_std(2,2:4),'ro','LineWidth',3)
errorbar(2:4,rms_ratio_mean(3,2:4),rms_ratio_std(3,2:4),'gx','LineWidth',3)
set(gca,'FontSize',16)
% axis tight
grid on
xticks([2:4])

%%  Monte Carlo plots for GPM 
% cases to denote for the paper: {A0} (Case 1), {A0, A4} (Case 2),
% {A0,A4,A8} (Case 3), {A0,A4,A8,A12} (Case 5), {A0,A8} (Case 4)
str_case = 'MC_bff_case';
ind_case = [1,2,3,5,4]; %[1,2,3,4,5,6,7,8]; %
for ii = 1:numel(ind_case)
    load(strcat(str_case,num2str(ii),'_mod'))
    cd_mean(1,ind_case(ii)) = mean(mc_cderr);
    cd_std(1,ind_case(ii)) = std(mc_cderr);
    pos_mean(1,ind_case(ii)) = mean(mc_poserr);
    pos_std(1,ind_case(ii)) = std(mc_poserr);   
    vel_mean(1,ind_case(ii)) = mean(mc_velerr);
    vel_std(1,ind_case(ii)) = std(mc_velerr);   
    rms_ratio_mean(1,ind_case(ii)) = mean(mc_rms_res_norm);
    rms_ratio_std(1,ind_case(ii)) = std(mc_rms_res_norm);
    
    load(strcat(str_case,num2str(ii)))
    cd_mean(2,ind_case(ii)) = mean(mc_cderr);
    cd_std(2,ind_case(ii)) = std(mc_cderr);
    pos_mean(2,ind_case(ii)) = mean(mc_poserr);
    pos_std(2,ind_case(ii)) = std(mc_poserr);   
    vel_mean(2,ind_case(ii)) = mean(mc_velerr);
    vel_std(2,ind_case(ii)) = std(mc_velerr);   
    rms_ratio_mean(2,ind_case(ii)) = mean(mc_rms_res_norm);
    rms_ratio_std(2,ind_case(ii)) = std(mc_rms_res_norm); 
    
    load(strcat(str_case,num2str(ii), '_mod_noErrOFF'))
    cd_mean(3,ind_case(ii)) = mean(mc_cderr);
    cd_std(3,ind_case(ii)) = std(mc_cderr);
    pos_mean(3,ind_case(ii)) = mean(mc_poserr);
    pos_std(3,ind_case(ii)) = std(mc_poserr);   
    vel_mean(3,ind_case(ii)) = mean(mc_velerr);
    vel_std(3,ind_case(ii)) = std(mc_velerr);   
    rms_ratio_mean(3,ind_case(ii)) = mean(mc_rms_res_norm);
    rms_ratio_std(3,ind_case(ii)) = std(mc_rms_res_norm);      
end
 
figure(1)
errorbar(cd_mean(1,:),cd_std(1,:),'k*','LineWidth',4)
hold on
errorbar(cd_mean(2,:),cd_std(2,:),'ro','LineWidth',1)
errorbar(cd_mean(3,:),cd_std(3,:),'gx','LineWidth',3)
ylabel('Drag-coefficient error')
xlabel('Case')
% set(gca,'XTickLabel', {'Case 1','Case 2','Case 3','Case 4','Case 5','Case 6','Case 7'});
title('Drag-coefficient error')
legend('Modeled higher-orders','Ignored higher-orders','True OFF coefficients')
set(gca,'FontSize',18)
% set(gca,'YScale','log')
grid on

% axes('position',[.3 .4 .47 .25])
% box on 
% errorbar(2:8,cd_mean(1,2:8),cd_std(1,2:8),'k*','LineWidth',4)
% hold on
% errorbar(2:8,cd_mean(2,2:8),cd_std(2,2:8),'ro','LineWidth',1)
% set(gca,'FontSize',16)
% % axis tight
% grid on
% ylim([0.01 .2])
% xticks([2:8])

figure(2)
errorbar(pos_mean(1,:),pos_std(1,:),'k*','LineWidth',4)
hold on
errorbar(pos_mean(2,:),pos_std(2,:),'ro','LineWidth',1)
errorbar(pos_mean(3,:),pos_std(3,:),'gx','LineWidth',3)
ylabel('Initial position error (m)')
% set(gca,'XTickLabel', {'Case 1','Case 2','Case 3','Case 4','Case 5','Case 6','Case 7'});
title('Initial position error')
legend('Modeled higher-orders','Ignored higher-orders','True OFF coefficients')
set(gca,'FontSize',18)
set(gca,'YScale','log')
grid on

% axes('position',[.3 .4 .5 .25])
% box on 
% errorbar(2:8,pos_mean(1,2:8),pos_std(1,2:8),'k*','LineWidth',3)
% hold on
% errorbar(2:8,pos_mean(2,2:8),pos_std(2,2:8),'ro','LineWidth',3)
% set(gca,'FontSize',16)
% % axis tight
% grid on
% xticks([2:8])
% ylim([0 0.6])

figure(3)
errorbar(rms_ratio_mean(1,:),rms_ratio_std(1,:),'k*','LineWidth',4)
hold on
errorbar(rms_ratio_mean(2,:),rms_ratio_std(2,:),'ro','LineWidth',1)
errorbar(rms_ratio_mean(3,:),rms_ratio_std(3,:),'gx','LineWidth',2)
ylabel('Measurement residual ratio')
% set(gca,'XTickLabel', {'Case 1','Case 2','Case 3','Case 4','Case 5','Case 6','Case 7'});
title('Measurement residual ratio')
legend('Modeled higher-orders','Ignored higher-orders','True OFF coefficients')
set(gca,'FontSize',18)
set(gca,'YScale','log')
grid on
% xlim([2 8])

axes('position',[.3 .4 .5 .25])
box on 
errorbar(2:4,rms_ratio_mean(1,2:4),rms_ratio_std(1,2:4),'k*','LineWidth',4)
hold on
errorbar(2:4,rms_ratio_mean(2,2:4),rms_ratio_std(2,2:4),'ro','LineWidth',1)
errorbar(2:4,rms_ratio_mean(3,2:4),rms_ratio_std(3,2:4),'gx','LineWidth',3)
set(gca,'FontSize',16)
set(gca,'YScale','log')
% axis tight
grid on
ylim([0.99 1.01])
xticks([2:4])