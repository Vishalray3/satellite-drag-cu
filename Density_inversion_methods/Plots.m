%% Plots
load Cd_sphere_dria_150F10
% Cd variation
figure(1)
plot(Cdtotal_mat(1,:), 'k','LineWidth', 1)
hold on
plot(Cdtotal_mat(end,:), '--r','LineWidth', 1)
xlabel('Eccentric anomaly (deg)')
ylabel('Drag-coefficient')
legend('K = 1e6', 'K = 1e8')
title('Drag-coefficient in orbit')
set(gca,'FontSize',18)
%%
figure(2)
ax1 = subplot(2,1,1);
plot(Kl_mat, Af_total_mat(1,:), 'k','LineWidth', 1)
title('Variation of Fourier coefficients')
leg = legend('$\bar{A}_0$');
set(leg,'interpreter','latex');
set(gca,'FontSize',18)
ax2 = subplot(2,1,2);
plot(Kl_mat, Af_total_mat(2,:), 'b','LineWidth', 1)
hold on
plot(Kl_mat, Af_total_mat(3,:), '--r','LineWidth', 1)
plot(Kl_mat, Af_total_mat(4,:), '-.c','LineWidth', 1)
plot(Kl_mat, Af_total_mat(5,:), ':m','LineWidth', 1)
leg = legend('$\bar{A}_1$','$\bar{A}_2$','$\bar{A}_3$','$\bar{A}_4$');
set(leg,'interpreter','latex');
xlabel('Langmuir adsorbate constant')
set(gca,'FontSize',18)

p1=get(ax1,'position');
p2=get(ax2,'position');
height=p1(2)+p1(4)-p2(2);
h3=axes('position',[p2(1) p2(2) p2(3) height],'visible','off');
h_label=ylabel('Fourier coefficients','visible','on');
set(gca,'FontSize',18)

%%
load Case1_sphere
xvec = time_prop_utc/T_orb;

figure(1)
plot(xvec, obs_rank(1,:), 'k','LineWidth', 1)
xlabel('Orbits')
ylabel('Rank')
title('Rank of observability matrix')
grid on
set(gca,'FontSize',18)
set(gca,'XTick',xvec(1):2:xvec(end))

figure(2)
plot(xvec, s_obs(end,:), 'k','LineWidth', 1)
hold on
plot(xvec, s_obs(end-1,:), '--r','LineWidth', 1)
plot(xvec, s_thresh, '-.g','LineWidth', 1)
leg = legend('$s_1$','$s_2$','threshold');
set(leg,'interpreter','latex')
xlabel('Orbits')
grid on
ylabel('Singular values')
title('Singular values of observability matrix')
set(gca,'FontSize',18)
set(gca,'XTick',xvec(1):2:xvec(end))

figure(3)
ax1 = subplot(2,1,1);
semilogy(xvec, sigma_x(7,:)/abs(X_nom(7)), 'k','LineWidth', 1)
hold on
semilogy(xvec, sigma_x(8,:)/abs(X_nom(8)), '--r','LineWidth', 1)
title('Uncertainty of Fourier coefficients')
leg = legend('$\bar{A}_0$','$\bar{A}_1$');
set(leg,'interpreter','latex');
set(gca,'FontSize',18)
% ylim([0 0.2])
grid on
ax2 = subplot(2,1,2);
plot(xvec, sigma_x(9,:)/abs(X_nom(9)), 'k','LineWidth', 1)
hold on
plot(xvec, sigma_x(10,:)/abs(X_nom(10)), '--r','LineWidth', 1)
plot(xvec, sigma_x(11,:)/abs(X_nom(11)), '-.g','LineWidth', 1)
leg = legend('$\bar{A}_2$','$\bar{A}_3$','$\bar{A}_4$');
set(leg,'interpreter','latex');
xlabel('Orbits')
set(gca,'FontSize',18)
ylim([0 3])
grid on
p1=get(ax1,'position');
p2=get(ax2,'position');
height=p1(2)+p1(4)-p2(2);
h3=axes('position',[p2(1) p2(2) p2(3) height],'visible','off');
h_label=ylabel('Standard deviation (ratio)','visible','on');
set(gca,'FontSize',18)

%% consider covariance results

load Case1_sphere_consider_order0_cbar
sig_pos_est = [norm(sigma_x(1:3,end)); norm(sigma_xx(1:3,end))];
sig_vel_est = [norm(sigma_x(4:6,end)); norm(sigma_xx(4:6,end))];
sig_cd_est = [norm(sigma_x(7,end)); norm(sigma_xx(7,end))];
load Case1_sphere_consider_order1_cbar
sig_pos_est = [sig_pos_est,[norm(sigma_x(1:3,end)); norm(sigma_xx(1:3,end))] ];
sig_vel_est = [sig_vel_est,[norm(sigma_x(4:6,end)); norm(sigma_xx(4:6,end))]];
sig_cd_est = [sig_cd_est,[norm(sigma_x(7,end));norm(sigma_xx(7,end))]];
load Case1_sphere_consider_order2_cbar
sig_pos_est = [sig_pos_est,[norm(sigma_x(1:3,end)); norm(sigma_xx(1:3,end))] ];
sig_vel_est = [sig_vel_est,[norm(sigma_x(4:6,end)); norm(sigma_xx(4:6,end))]];
sig_cd_est = [sig_cd_est,[norm(sigma_x(7,end));norm(sigma_xx(7,end))]];
load Case1_sphere_consider_order3_cbar
sig_pos_est = [sig_pos_est,[norm(sigma_x(1:3,end)); norm(sigma_xx(1:3,end))] ];
sig_vel_est = [sig_vel_est,[norm(sigma_x(4:6,end)); norm(sigma_xx(4:6,end))]];
sig_cd_est = [sig_cd_est,[norm(sigma_x(7,end));norm(sigma_xx(7,end))]];
load Case1_sphere_consider_order4_cbar
sig_pos_est = [sig_pos_est,[norm(sigma_x(1:3,end)); norm(sigma_xx(1:3,end))] ];
sig_vel_est = [sig_vel_est,[norm(sigma_x(4:6,end)); norm(sigma_xx(4:6,end))]];
sig_cd_est = [sig_cd_est,[norm(sigma_x(7,end));norm(sigma_xx(7,end))]];

vec = 0:4;
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

%% Correlation matrix
figure()
imagesc(abs(Corr_mat(7:end,7:end)))
colormap copper
c = colorbar;
c.Label.String = 'Correlation coefficients';
title('Correlation matrix')
xticks([1 2 3 4 5 6 7 8 9 10 11])
yticks([1 2 3 4 5 6 7 8 9 10 11])
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTickLabel', {'$\overline{A}_0$','$\overline{A}_1$','$\overline{A}_2$','$\overline{A}_3$',...
    '$\overline{A}_4$','$\overline{A}_5$','$\overline{A}_6$','$\overline{A}_7$','$\overline{A}_8$',...
    '$\overline{A}_9$','$\overline{A}_{10}$'});
set(gca,'YTickLabel', {'$\overline{A}_0$','$\overline{A}_1$','$\overline{A}_2$','$\overline{A}_3$',...
    '$\overline{A}_4$','$\overline{A}_5$','$\overline{A}_6$','$\overline{A}_7$','$\overline{A}_8$',...
    '$\overline{A}_9$','$\overline{A}_{10}$'});
set(gca,'FontSize',18)

%% Plot consider covariance errors
perr = [];
cderr = [];
resp = [];
load Case1_sphere_order0
perr = [perr,pos_err];
cderr = [cderr,Cd_err_rms];
resp = [resp,norm(rms_res_post(1:3,end))];
load Case1_sphere_order1
perr = [perr,pos_err];
cderr = [cderr,Cd_err_rms];
resp = [resp,norm(rms_res_post(1:3,end))];
load Case1_sphere_order2
perr = [perr,pos_err];
cderr = [cderr,Cd_err_rms];
resp = [resp,norm(rms_res_post(1:3,end))];
load Case1_sphere_order3
perr = [perr,pos_err];
cderr = [cderr,Cd_err_rms];
resp = [resp,norm(rms_res_post(1:3,end))];
load Case1_sphere_order4
perr = [perr,pos_err];
cderr = [cderr,Cd_err_rms];
resp = [resp,norm(rms_res_post(1:3,end))];
% load cbartrue_sphere_consider_order5
% perr = [perr,pos_err];
% cderr = [cderr,Cd_err_rms];
% resp = [resp,norm(rms_res_post(1:3,end))];

err_mat = [perr;cderr;resp];
xlab = [{'Init. pos. (m)'},{'C_d RMS'},{'Pos. res. RMS (m)'}];
figure(3)
bar(err_mat)
set(gca,'yscale','log')
xticklabels(xlab)
set(gca,'FontSize',18)
% xtickangle(45)
ylabel('Error norm')
legend('Order 0','Order 1','Order 2','Order 3','Order 4','Order 5')
title('Truncation errors')

%%
[sig_rat, ind_rat] = sort(sigma_x(7:end,1)./sigma_x(7:end,end), 'descend');
ind_rat = ind_rat-1;
figure()
bar(sig_rat)
% set(gca,'yscale','log')
set(gca,'FontSize',18)
ylabel('Standard deviation ratio')
xlabel('Ranked Fourier coefficients')
title('Standard deviation ratio')
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTickLabel', {'$\overline{A}_0$','$\overline{A}_1$','$\overline{A}_2$','$\overline{A}_3$',...
    '$\overline{A}_{10}$','$\overline{A}_{4}$','$\overline{A}_5$','$\overline{A}_6$','$\overline{A}_7$',...
    '$\overline{A}_{9}$','$\overline{A}_{8}$'});
% ylim([1e-10 1e7])
grid on
%% Monte Carlo plots 
% cases to denote for the paper: {A0} (Case 6), {A0, A1} (Case 4),
% {A0,A1,A2} (Case 7), {A0,A1,A2,A3} (Case 1), {A0,A1,A2,A3,A4} (Case 5),
% {A0,A1,A3} (Case 3),{A0,A2,A3} (Case 2)
str_case = 'MC_case';
ind_case = [4,7,6,2,5,1,3];
for ii = 1:7
    load(strcat(str_case,num2str(ii)))
    cd_mean(1,ind_case(ii)) = mean(mc_cderr);
    cd_std(1,ind_case(ii)) = std(mc_cderr);
    pos_mean(1,ind_case(ii)) = mean(mc_poserr);
    pos_std(1,ind_case(ii)) = std(mc_poserr);   
    vel_mean(1,ind_case(ii)) = mean(mc_velerr);
    vel_std(1,ind_case(ii)) = std(mc_velerr);   
    rms_ratio_mean(1,ind_case(ii)) = mean(mc_rms_res_norm);
    rms_ratio_std(1,ind_case(ii)) = std(mc_rms_res_norm);
    
    load(strcat(str_case,num2str(ii), '_unmod'))
    cd_mean(2,ind_case(ii)) = mean(mc_cderr);
    cd_std(2,ind_case(ii)) = std(mc_cderr);
    pos_mean(2,ind_case(ii)) = mean(mc_poserr);
    pos_std(2,ind_case(ii)) = std(mc_poserr);   
    vel_mean(2,ind_case(ii)) = mean(mc_velerr);
    vel_std(2,ind_case(ii)) = std(mc_velerr);   
    rms_ratio_mean(2,ind_case(ii)) = mean(mc_rms_res_norm);
    rms_ratio_std(2,ind_case(ii)) = std(mc_rms_res_norm);    
end
 
figure(1)
errorbar(cd_mean(1,:),cd_std(1,:),'k*','LineWidth',4)
hold on
errorbar(cd_mean(2,:),cd_std(2,:),'ro','LineWidth',1)
ylabel('Drag-coefficient error')
xlabel('Case')
% set(gca,'XTickLabel', {'Case 1','Case 2','Case 3','Case 4','Case 5','Case 6','Case 7'});
title('Drag-coefficient error')
legend('Modeled higher-orders','Ignored higher-orders')
set(gca,'FontSize',18)
set(gca,'YScale','log')
grid on

% axes('position',[.3 .4 .47 .25])
% box on 
% errorbar(2:6,cd_mean(1,2:6),cd_std(1,2:6),'k*','LineWidth',1)
% hold on
% errorbar(2:6,cd_mean(2,2:6),cd_std(2,2:6),'ro','LineWidth',1)
% set(gca,'FontSize',16)
% % axis tight
% grid on
% 
figure(2)
errorbar(pos_mean(1,:),pos_std(1,:),'k*','LineWidth',4)
hold on
errorbar(pos_mean(2,:),pos_std(2,:),'ro','LineWidth',1)
ylabel('Initial position error (m)')
% set(gca,'XTickLabel', {'Case 1','Case 2','Case 3','Case 4','Case 5','Case 6','Case 7'});
title('Initial position error')
legend('Modeled higher-orders','Ignored higher-orders')
set(gca,'FontSize',18)
set(gca,'YScale','log')
grid on
% 
% axes('position',[.3 .4 .5 .25])
% box on 
% errorbar(2:7,pos_mean(1,2:7),pos_std(1,2:7),'k*','LineWidth',1)
% hold on
% errorbar(2:7,pos_mean(2,2:7),pos_std(2,2:7),'ro','LineWidth',1)
% set(gca,'FontSize',16)
% % axis tight
% grid on
% xticks([2:7])
% 
% figure(3)
% errorbar(rms_ratio_mean(1,:),rms_ratio_std(1,:),'k*','LineWidth',1)
% hold on
% errorbar(rms_ratio_mean(2,:),rms_ratio_std(2,:),'ro','LineWidth',1)
% ylabel('Drag-coefficient error')
% % set(gca,'XTickLabel', {'Case 1','Case 2','Case 3','Case 4','Case 5','Case 6','Case 7'});
% title('Monte Carlo results')
% legend('Modeled higher-orders','Ignored higher-orders')
% set(gca,'FontSize',18)
% grid on