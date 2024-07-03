%% PLots bodf
%% Standard deviation
%GPM
% xvec = time_prop_utc/T_orb;
% figure(3)
% ax1 = subplot(2,1,1);
% plot(xvec, sigma_x(7,:)/abs(X_nom(7)), 'k','LineWidth', 1)
% hold on
% semilogy(xvec, sigma_x(8,:)/abs(X_nom(8)), 'r','LineWidth', 1)
% semilogy(xvec, sigma_x(15,:)/abs(X_nom(15)), 'b','LineWidth', 1)
% semilogy(xvec, sigma_x(26,:)/abs(X_nom(26)), 'c','LineWidth', 1)
% semilogy(xvec, sigma_x(23,:)/abs(X_nom(23)), 'g','LineWidth', 1)
% semilogy(xvec, sigma_x(31,:)/abs(X_nom(31)), 'm','LineWidth', 1)
% title('Uncertainty of Fourier coefficients')
% leg = legend('$\overline{A}_{0,0}$','$\overline{A}_{2,0}$','$\overline{A}_{4,0}$','$\overline{A}_{0,1}$','$\overline{B}_{1,0}$');
% set(leg,'interpreter','latex');
% set(gca,'FontSize',18)
% % set(gca,'YTick',[1e-3,1e-2,1e-1])
% grid on
% xlim([0 17])
% ax2 = subplot(2,1,2);
% plot(xvec, sigma_x(13,:)/abs(X_nom(13)), 'k','LineWidth', 1)
% hold on
% semilogy(xvec, sigma_x(30,:)/abs(X_nom(30)), 'r','LineWidth', 1)
% semilogy(xvec, sigma_x(68,:)/abs(X_nom(68)), 'b','LineWidth', 1)
% semilogy(xvec, sigma_x(97,:)/abs(X_nom(97)), 'c','LineWidth', 1)
% leg = legend('$\overline{A}_{6,0}$','$\overline{A}_{2,1}$','$\overline{A}_{2,3}$','$\overline{B}_{1,1}$');
% set(leg,'interpreter','latex');
% xlabel('Orbits')
% set(gca,'FontSize',18)
% % ylim([0 0.2])
% grid on
% xlim([0 17])
% p1=get(ax1,'position');
% p2=get(ax2,'position');
% height=p1(2)+p1(4)-p2(2);
% h3=axes('position',[p2(1) p2(2) p2(3) height],'visible','off');
% h_label=ylabel('Standard deviation','visible','on');
% set(gca,'FontSize',18)

%cube
xvec = time_prop_utc/T_orb;
figure(3)
ax1 = subplot(2,1,1);
plot(xvec, sigma_x(7,:)/abs(X_nom(7)), 'k','LineWidth', 1)
hold on
semilogy(xvec, sigma_x(8,:)/abs(X_nom(8)), 'r','LineWidth', 1)
semilogy(xvec, sigma_x(15,:)/abs(X_nom(15)), 'b','LineWidth', 1)
title('Uncertainty of Fourier coefficients')
leg = legend('$\overline{A}_{0,0}$','$\overline{A}_{4,0}$','$\overline{A}_{0,1}$');
set(leg,'interpreter','latex');
set(gca,'FontSize',18)
% set(gca,'YTick',[1e-3,1e-2,1e-1])
grid on
xlim([0 17])
ax2 = subplot(2,1,2);
semilogy(xvec, sigma_x(26,:)/abs(X_nom(26)), 'k','LineWidth', 1)
hold on
semilogy(xvec, sigma_x(23,:)/abs(X_nom(23)), 'r','LineWidth', 1)
semilogy(xvec, sigma_x(31,:)/abs(X_nom(31)), 'b','LineWidth', 1)
leg = legend('$\overline{A}_{4,1}$','$\overline{A}_{0,2}$','$\overline{A}_{0,3}$');
set(leg,'interpreter','latex');
xlabel('Orbits')
set(gca,'FontSize',18)
% ylim([0 0.2])
grid on
xlim([0 17])
p1=get(ax1,'position');
p2=get(ax2,'position');
height=p1(2)+p1(4)-0.2*p2(2);
h3=axes('position',[p2(1) p2(2) p2(3) height],'visible','off');
h_label=ylabel('Standard deviation','visible','on');
set(gca,'FontSize',18)

%% Correlation matrix
% gpm
% figure()
% imagesc(abs(Corr_mat([7,9,11,13,28,30,68,92,97],[7,9,11,13,28,30,68,92,97])))
% colormap copper
% c = colorbar;
% c.Label.String = 'Correlation coefficients';
% title('Correlation matrix')
% xticks([1:9])
% yticks([1:9])
% set(gca,'TickLabelInterpreter', 'latex');
% set(gca,'XTickLabel', {'$\overline{A}_{0,0}$','$\overline{A}_{2,0}$','$\overline{A}_{4,0}$','$\overline{A}_{6,0}$','$\overline{A}_{0,1}$','$\overline{A}_{2,1}$',...
%    '$\overline{A}_{2,3}$','$\overline{B}_{1,0}$','$\overline{B}_{1,1}$'});
% set(gca,'YTickLabel', {'$\overline{A}_{0,0}$','$\overline{A}_{2,0}$','$\overline{A}_{4,0}$','$\overline{A}_{6,0}$','$\overline{A}_{0,1}$','$\overline{A}_{2,1}$',...
%    '$\overline{A}_{2,3}$','$\overline{B}_{1,0}$','$\overline{B}_{1,1}$'});
% set(gca,'FontSize',18)

% cube
figure()
imagesc(abs(Corr_mat([7,8,15,23,26,31],[7,8,15,23,26,31])))
colormap copper
c = colorbar;
c.Label.String = 'Correlation coefficients';
title('Correlation matrix')
xticks([1:6])
yticks([1:6])
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTickLabel', {'$\overline{A}_{0,0}$','$\overline{A}_{4,0}$','$\overline{A}_{0,1}$','$\overline{A}_{4,1}$','$\overline{A}_{0,2}$','$\overline{A}_{0,3}$'});
set(gca,'YTickLabel', {'$\overline{A}_{0,0}$','$\overline{A}_{4,0}$','$\overline{A}_{0,1}$','$\overline{A}_{4,1}$','$\overline{A}_{0,2}$','$\overline{A}_{0,3}$'});
set(gca,'FontSize',18)
%% Consider covariance
load bodf_cube_consider_case1
sig_pos_est = [norm(sigma_x(1:3,end)); norm(sigma_xx(1:3,end))];
sig_vel_est = [norm(sigma_x(4:6,end)); norm(sigma_xx(4:6,end))];
sig_cd_est = [norm(sigma_x(7,end)); norm(sigma_xx(7,end))];
load bodf_cube_consider_case2
sig_pos_est = [sig_pos_est,[norm(sigma_x(1:3,end)); norm(sigma_xx(1:3,end))] ];
sig_vel_est = [sig_vel_est,[norm(sigma_x(4:6,end)); norm(sigma_xx(4:6,end))]];
sig_cd_est = [sig_cd_est,[norm(sigma_x(7,end));norm(sigma_xx(7,end))]];
load bodf_cube_consider_case3
sig_pos_est = [sig_pos_est,[norm(sigma_x(1:3,end)); norm(sigma_xx(1:3,end))] ];
sig_vel_est = [sig_vel_est,[norm(sigma_x(4:6,end)); norm(sigma_xx(4:6,end))]];
sig_cd_est = [sig_cd_est,[norm(sigma_x(7,end));norm(sigma_xx(7,end))]];
load bodf_cube_consider_case4
sig_pos_est = [sig_pos_est,[norm(sigma_x(1:3,end)); norm(sigma_xx(1:3,end))] ];
sig_vel_est = [sig_vel_est,[norm(sigma_x(4:6,end)); norm(sigma_xx(4:6,end))]];
sig_cd_est = [sig_cd_est,[norm(sigma_x(7,end));norm(sigma_xx(7,end))]];
load bodf_cube_consider_case5
sig_pos_est = [sig_pos_est,[norm(sigma_x(1:3,end)); norm(sigma_xx(1:3,end))] ];
sig_vel_est = [sig_vel_est,[norm(sigma_x(4:6,end)); norm(sigma_xx(4:6,end))]];
sig_cd_est = [sig_cd_est,[norm(sigma_x(7,end));norm(sigma_xx(7,end))]];
load bodf_cube_consider_case6
sig_pos_est = [sig_pos_est,[norm(sigma_x(1:3,end)); norm(sigma_xx(1:3,end))] ];
sig_vel_est = [sig_vel_est,[norm(sigma_x(4:6,end)); norm(sigma_xx(4:6,end))]];
sig_cd_est = [sig_cd_est,[norm(sigma_x(7,end));norm(sigma_xx(7,end))]];
% load bodf_gpm_consider_case7
% sig_pos_est = [sig_pos_est,[norm(sigma_x(1:3,end)); norm(sigma_xx(1:3,end))] ];
% sig_vel_est = [sig_vel_est,[norm(sigma_x(4:6,end)); norm(sigma_xx(4:6,end))]];
% sig_cd_est = [sig_cd_est,[norm(sigma_x(7,end));norm(sigma_xx(7,end))]];
% load bodf_gpm_consider_case8
% sig_pos_est = [sig_pos_est,[norm(sigma_x(1:3,end)); norm(sigma_xx(1:3,end))] ];
% sig_vel_est = [sig_vel_est,[norm(sigma_x(4:6,end)); norm(sigma_xx(4:6,end))]];
% sig_cd_est = [sig_cd_est,[norm(sigma_x(7,end));norm(sigma_xx(7,end))]];
vec = [1:6];
figure()
semilogy(vec, sig_pos_est(1,:),'-k','LineWidth',1)
hold on
semilogy(vec, sig_vel_est(1,:),'-r','LineWidth',1)
semilogy(vec, sig_cd_est(1,:),'-g','LineWidth',1)
semilogy(vec, sig_pos_est(2,:),'--k','LineWidth',1)
semilogy(vec, sig_vel_est(2,:),'--r','LineWidth',1)
semilogy(vec, sig_cd_est(2,:),'--g','LineWidth',1)
legend({'Position','Velocity','Order 0'})
xlabel('Cases')
ylabel('Standard deviation')
title('Estimated and considered uncertainties')
set(gca,'FontSize',18)
set(gca,'XTick',vec)

%%  Monte Carlo plots for GPM 
% cases to denote for the paper: {A0} (Case 1), {A0, A4} (Case 2),
% {A0,A4,A8} (Case 3), {A0,A4,A8,A12} (Case 5), {A0,A8} (Case 4)
str_case = 'MC_bodf_case';
ind_case_mc = [1:9]; %[1,2,3,4,5,6]; 
for ii = 1:numel(ind_case_mc)
    load(strcat(str_case,num2str(ii),'_gpm_mod'))
    cd_mean(1,ind_case_mc(ii)) = mean(mc_cderr);
    cd_std(1,ind_case_mc(ii)) = std(mc_cderr);
    pos_mean(1,ind_case_mc(ii)) = mean(mc_poserr);
    pos_std(1,ind_case_mc(ii)) = std(mc_poserr);   
    vel_mean(1,ind_case_mc(ii)) = mean(mc_velerr);
    vel_std(1,ind_case_mc(ii)) = std(mc_velerr);   
    rms_ratio_mean(1,ind_case_mc(ii)) = mean(mc_rms_res_norm);
    rms_ratio_std(1,ind_case_mc(ii)) = std(mc_rms_res_norm);
    
    load(strcat(str_case,num2str(ii),'_gpm'))
    cd_mean(2,ind_case_mc(ii)) = mean(mc_cderr);
    cd_std(2,ind_case_mc(ii)) = std(mc_cderr);
    pos_mean(2,ind_case_mc(ii)) = mean(mc_poserr);
    pos_std(2,ind_case_mc(ii)) = std(mc_poserr);   
    vel_mean(2,ind_case_mc(ii)) = mean(mc_velerr);
    vel_std(2,ind_case_mc(ii)) = std(mc_velerr);   
    rms_ratio_mean(2,ind_case_mc(ii)) = mean(mc_rms_res_norm);
    rms_ratio_std(2,ind_case_mc(ii)) = std(mc_rms_res_norm); 
    
%     load(strcat(str_case,num2str(ii), '_mod_noErrOFF'))
%     cd_mean(3,ind_case(ii)) = mean(mc_cderr);
%     cd_std(3,ind_case(ii)) = std(mc_cderr);
%     pos_mean(3,ind_case(ii)) = mean(mc_poserr);
%     pos_std(3,ind_case(ii)) = std(mc_poserr);   
%     vel_mean(3,ind_case(ii)) = mean(mc_velerr);
%     vel_std(3,ind_case(ii)) = std(mc_velerr);   
%     rms_ratio_mean(3,ind_case(ii)) = mean(mc_rms_res_norm);
%     rms_ratio_std(3,ind_case(ii)) = std(mc_rms_res_norm);      
end
 
vec_box = 2:6;
figure(1)
errorbar(cd_mean(1,:),cd_std(1,:),'k*','LineWidth',4)
hold on
errorbar(cd_mean(2,:),cd_std(2,:),'ro','LineWidth',1)
% errorbar(cd_mean(3,:),cd_std(3,:),'gx','LineWidth',3)
ylabel('Drag-coefficient error')
xlabel('Case')
% set(gca,'XTickLabel', {'Case 1','Case 2','Case 3','Case 4','Case 5','Case 6','Case 7'});
title('Drag-coefficient error')
legend('Modeled higher-orders','Ignored higher-orders','True OFF coefficients')
set(gca,'FontSize',18)
set(gca,'YScale','log')
grid on

% axes('position',[.3 .4 .47 .25])
% box on 
% errorbar(vec_box,cd_mean(1,vec_box),cd_std(1,vec_box),'k*','LineWidth',3)
% hold on
% errorbar(vec_box,cd_mean(2,vec_box),cd_std(2,vec_box),'ro','LineWidth',3)
% set(gca,'FontSize',16)
% % axis tight
% grid on
% ylim([0.0 .2])
% xticks(vec_box)

figure(2)
errorbar(pos_mean(1,:),pos_std(1,:),'k*','LineWidth',4)
hold on
errorbar(pos_mean(2,:),pos_std(2,:),'ro','LineWidth',1)
% errorbar(pos_mean(3,:),pos_std(3,:),'gx','LineWidth',3)
ylabel('Initial position error (m)')
xlabel('Case')
% set(gca,'XTickLabel', {'Case 1','Case 2','Case 3','Case 4','Case 5','Case 6','Case 7'});
title('Initial position error')
legend('Modeled higher-orders','Ignored higher-orders','True OFF coefficients')
set(gca,'FontSize',18)
set(gca,'YScale','log')
grid on

% axes('position',[.3 .4 .5 .25])
% box on 
% errorbar(vec_box,pos_mean(1,vec_box),pos_std(1,vec_box),'k*','LineWidth',3)
% hold on
% errorbar(vec_box,pos_mean(2,vec_box),pos_std(2,vec_box),'ro','LineWidth',3)
% set(gca,'FontSize',16)
% % axis tight
% grid on
% xticks(vec_box)
% % ylim([0 0.6])

% vec_box = 4:9;
figure(3)
errorbar(rms_ratio_mean(1,:),rms_ratio_std(1,:),'k*','LineWidth',4)
hold on
errorbar(rms_ratio_mean(2,:),rms_ratio_std(2,:),'ro','LineWidth',1)
% errorbar(rms_ratio_mean(3,:),rms_ratio_std(3,:),'gx','LineWidth',3)
ylabel('Measurement residual ratio')
xlabel('Case')
% set(gca,'XTickLabel', {'Case 1','Case 2','Case 3','Case 4','Case 5','Case 6','Case 7'});
title('Measurement residual ratio')
legend('Modeled higher-orders','Ignored higher-orders','True OFF coefficients')
set(gca,'FontSize',18)
set(gca,'YScale','log')
grid on

% vec_box = [4:9]
% axes('position',[.3 .4 .5 .25])
% box on 
% errorbar(vec_box,rms_ratio_mean(1,vec_box),rms_ratio_std(1,vec_box),'k*','LineWidth',4)
% hold on
% errorbar(vec_box,rms_ratio_mean(2,vec_box),rms_ratio_std(2,vec_box),'ro','LineWidth',1)
% % errorbar(2:4,rms_ratio_mean(3,2:4),rms_ratio_std(3,2:4),'gx','LineWidth',3)
% set(gca,'FontSize',16)
% % axis tight
% grid on
% % ylim([0.9998 1.0008])
% xticks(vec_box)
