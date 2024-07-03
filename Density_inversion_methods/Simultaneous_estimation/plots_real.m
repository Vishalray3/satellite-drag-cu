%% plotting real data
load 1_ekf_msis_cr_cdorder2
rho_est1 = rho_est;
Cds_est1 = Cds_est;
Xs_est1 = Xs_est;

load 1_ekf_jb08_cr_cdorder2

%%
close all
clc
clearvars
% figure(2)
% plot(time_prop_utc/T_orb, rho_est1 - rho_est)
% hold on
% plot(time_prop_utc/T_orb,rho_est.*(Xs_est(10,:))-rho_est1.*(Xs_est1(10,:)))
% plot(time_prop_utc/T_orb,rho_est.*(Xs_est(10,:) + Xs_est(8,:))-rho_est1.*(Xs_est1(10,:)+ Xs_est1(8,:)))
load 2_ekf_msis_cr_cdorder2_zeta8_1e-7
figure(1)
plot(time_prop_utc/T_orb, rho_est, 'k','LineWidth', 1)
hold on
plot(time_prop_utc/T_orb,rho_est.*(Xs_est(7,:) + Xs_est(9,:))+rho_est, '--r','LineWidth', 1)
plot(time_prop_utc/T_orb,rho_true,'--b','LineWidth', 1)
xlim([2 15])
xlabel('Orbits')
ylabel('Density')
legend('NRLMSISE-00', 'Filter','HASDM')
title('Density estimate')
set(gca,'FontSize',18)
grid on

ind_vec= (2000:16000);
mean_rho0 = mean(rho_true(ind_vec) - rho_est(ind_vec))
mean_rho1 = mean(rho_true(ind_vec) - rho_est(ind_vec).*(Xs_est(7,(ind_vec)) + Xs_est(9,(ind_vec))))

rms_rho0 = rms(rho_true(ind_vec) - rho_est(ind_vec))
rms_rho1 = rms(rho_true(ind_vec) - rho_est(ind_vec).*(Xs_est(7,(ind_vec)) + Xs_est(9,(ind_vec))))

load 2_ekf_jb08_cr_cdorder2_zeta8_1e-7
figure(2)
plot(time_prop_utc/T_orb, rho_est, 'k','LineWidth', 1)
hold on
plot(time_prop_utc/T_orb,rho_est.*(Xs_est(7,:) + Xs_est(9,:))+rho_est, '--r','LineWidth', 1)
plot(time_prop_utc/T_orb,rho_true,'--b','LineWidth', 1)
xlim([2 15])
xlabel('Orbits')
ylabel('Density')
legend('JB2008', 'Filter','HASDM')
title('Density estimate')
set(gca,'FontSize',18)
grid on

mean_rho0 = mean(rho_true(ind_vec) - rho_est(ind_vec))
mean_rho1 = mean(rho_true(ind_vec) - rho_est(ind_vec).*(Xs_est(7,(ind_vec)) + Xs_est(9,(ind_vec))))

rms_rho0 = rms(rho_true(ind_vec) - rho_est(ind_vec))
rms_rho1 = rms(rho_true(ind_vec) - rho_est(ind_vec).*(Xs_est(7,(ind_vec)) + Xs_est(9,(ind_vec))))
% figure(1)
% plot(time_prop_utc/T_orb,rho_est.*(Xs_est(10,:) + Xs_est(8,:))+rho_est)
% hold on
% plot(time_prop_utc/T_orb,rho_est1.*(Xs_est1(10,:) + Xs_est1(8,:))+rho_est1)


%% Plots for propagation

figure(1)
subplot(2,1,1)
plot(sod_sec_sp3/3600, vecnorm(err_x_grav80_sp3(1:3,:),2,1),'.', 'MarkerSize',10)
hold on
plot(sod_sec_sp3/3600, vecnorm(err_x_third_sp3(1:3,:),2,1),'.', 'MarkerSize',10)
plot(sod_sec_sp3/3600, vecnorm(err_x_drag_msis_cball_sp3(1:3,:),2,1),'.', 'MarkerSize',10)
plot(sod_sec_sp3/3600, vecnorm(err_x_drag_jb08_cball_sp3(1:3,:),2,1),'.', 'MarkerSize',10)
plot(sod_sec_sp3/3600, vecnorm(err_x_drag_hasdm_cball_sp3(1:3,:),2,1),'.', 'MarkerSize',10)
ylabel('Error norm (m)')
title('Position error w.r.t POD')
leg = legend('Order 80 gravity','Third-body','Drag (MSIS00)','Drag (JB08)','Drag (HASDM)');
set(leg, 'interpreter','latex')
set(gca,'FontSize',14)
grid on

subplot(2,1,2)
plot(sod_sec_sp3/3600, vecnorm(err_x_drag_hasdm_cball_sp3(1:3,:),2,1),'.', 'MarkerSize',10)
hold on
plot(sod_sec_sp3/3600, vecnorm(err_x_drag_hasdmcball_tides_sp3(1:3,:),2,1),'.', 'MarkerSize',10)
plot(sod_sec_sp3/3600, vecnorm(err_x_drag_hasdmcball_tidessrp_sp3(1:3,:),2,1),'.', 'MarkerSize',10)
% plot(sod_sec_sp3/3600, vecnorm(err_x_drag_hasdmcball_all_sp3(1:3,:),2,1),'.', 'MarkerSize',10)
% plot(sod_sec_sp3/3600, vecnorm(err_x_drag_hasdmcball_120grav_all_sp3(1:3,:),2,1),'.', 'MarkerSize',10)
xlabel('Time (hours)')
ylabel('Error norm (m)')
leg = legend('Drag (HASDM)','Tides','SRP (Panel)');
set(leg, 'interpreter','latex')
set(gca,'FontSize',14)
grid on

%%
load spire83_propagation_2018_11_7
figure(1)
% subplot(2,1,1)
plot(sod_sec_sp3/3600, vecnorm(err_x_drag_hasdmcball_all_sp3(1:3,:),2,1),'.', 'MarkerSize',10)
hold on
plot(sod_sec_sp3/3600, vecnorm(err_x_drag_hasdm_paneldiffsm200_sp3(1:3,:),2,1),'.', 'MarkerSize',10)
plot(sod_sec_sp3/3600, vecnorm(err_x_drag_hasdmpanela93_all_sp3(1:3,:),2,1),'.', 'MarkerSize',10)
plot(sod_sec_sp3/3600, vecnorm(err_x_drag_hasdmpanelspec_all_sp3(1:3,:),2,1),'.', 'MarkerSize',10)
ylabel('Error norm (m)')
xlabel('Time (hours)')
title('Position error w.r.t POD (HASDM)')
leg = legend('Cannonball','DRIA, $f=1, \alpha = 1$', 'DRIA, $f=1, \alpha = 0.93$', 'DRIA, $f=0, \alpha = \alpha_s$');
set(leg, 'interpreter','latex', 'FontSize', 15)
set(gca,'FontSize',14)
grid on
ylim([0 1200])

figure(2)
% subplot(2,1,1)
plot(sod_sec_sp3/3600, vecnorm(err_x_drag_msis_cball_sp3(1:3,:),2,1),'.', 'MarkerSize',10)
hold on
plot(sod_sec_sp3/3600, vecnorm(err_x_drag_msispaneldiff_all_sp3(1:3,:),2,1),'.', 'MarkerSize',10)
plot(sod_sec_sp3/3600, vecnorm(err_x_drag_msispanela93_all_sp3(1:3,:),2,1),'.', 'MarkerSize',10)
plot(sod_sec_sp3/3600, vecnorm(err_x_drag_msispanelspec_all_sp3(1:3,:),2,1),'.', 'MarkerSize',10)
ylabel('Error norm (m)')
xlabel('Time (hours)')
title('Position error w.r.t POD (NRLMSISE-00)')
leg = legend('Cannonball','DRIA, $f=1, \alpha = 1$', 'DRIA, $f=1, \alpha = 0.93$', 'DRIA, $f=0, \alpha = \alpha_s$');
set(leg, 'interpreter','latex','FontSize', 15)
set(gca,'FontSize',14)
grid on
%%
load spire83_propagation_2018_11_7
grav80 = rms(vecnorm(err_x_grav80_sp3(1:3,:),2,1));
third = rms(vecnorm(err_x_third_sp3(1:3,:),2,1));
msis = rms(vecnorm(err_x_drag_msis_cball_sp3(1:3,:),2,1));
jb08 = rms(vecnorm(err_x_drag_jb08_cball_sp3(1:3,:),2,1));
hasdm = rms(vecnorm(err_x_drag_hasdm_cball_sp3(1:3,:),2,1));
tides = rms(vecnorm(err_x_drag_hasdmcball_tides_sp3(1:3,:),2,1));
srp = rms(vecnorm(err_x_drag_hasdmcball_tidessrp_sp3(1:3,:),2,1));
all = rms(vecnorm(err_x_drag_hasdmcball_all_sp3(1:3,:),2,1));
grav120 = rms(vecnorm(err_x_drag_hasdmcball_120grav_all_sp3(1:3,:),2,1));

err_prop = [grav80;third;msis;jb08;hasdm;tides;srp;all; grav120];

xlab =  [{'Order 80 gravity'},{'Third-body'},{'Drag (MSIS00)'},{'Drag (JB2008)'},{'Drag (HASDM)'}, {'Tides'}, {'SRP (Panel)'},{'Relativity'},{'Order 120 gravity'}];
figure()
bar(err_prop)
xticklabels(xlab)
xtickangle(30)
ylabel('RMS (m)')
xlabel('Force models')
title('RMS propagation errors')
set(gca,'FontSize',14)


%% Density and cd plots
load spire83_propagation_2018_11_7
figure(1)
% subplot(2,1,1)
plot(sod_uni(1:del_T:end)/3600, rho_msis, 'LineWidth',2)
hold on
plot(sod_uni(1:del_T:end)/3600, rho_jb08,':', 'LineWidth',2)
plot(sod_uni(1:del_T:end)/3600, rho_hasdm,'--k', 'LineWidth',2)
ylabel('Density ($kg/m^3$)','interpreter','latex')
xlabel('Time (hours)')
title('Density models')
leg = legend('NRLMSISE-00','JB2008','HASDM');
set(leg, 'interpreter','latex','FontSize', 15)
set(gca,'FontSize',14)
grid on

%%
load spire83_propagation_2018_11_7
load spire83_attitude_2018_11_7
vec = 1:numel(sod_uni);
figure(2)
% subplot(2,1,1)
plot(sod_uni(vec)/3600, Cd_diff(vec), 'LineWidth',1)
hold on
plot(sod_uni(vec)/3600, Cd_a93(vec) ,':', 'LineWidth',1)
plot(sod_uni(vec)/3600, Cd_spec(vec),'k--', 'LineWidth',1)
ylabel('Drag-coefficient')
xlabel('Time (hours)')
title('Drag-coefficient models')
leg = legend('DRIA, $f=1, \alpha = 1$', 'DRIA, $f=1, \alpha = 0.93$', 'DRIA, $f=0, \alpha = \alpha_s$');
set(leg, 'interpreter','latex','FontSize', 15)
set(gca,'FontSize',14)
grid on

%% Batch OD plots
load('spire83_batch_cdcrcball', 'y_res','sod_sec_sp3')
y_res1 = y_res;
load('spire83_batch_cdcballcrpanel', 'y_res')
figure(2)
subplot(2,1,1)
plot(sod_sec_sp3/3600, vecnorm(y_res1(1:3,:),2,1),'.', 'MarkerSize',10)
hold on
plot(sod_sec_sp3/3600, vecnorm(y_res(1:3,:),2,1),'.', 'MarkerSize',10)
ylabel('Position norm (m)')
title('Measurement residuals')
leg = legend('Cannonball SRP','Panel SRP');
set(leg, 'interpreter','latex')
set(gca,'FontSize',14)
grid on

subplot(2,1,2)
plot(sod_sec_sp3/3600, vecnorm(y_res1(4:6,:),2,1),'.', 'MarkerSize',10)
hold on
plot(sod_sec_sp3/3600, vecnorm(y_res(4:6,:),2,1),'.', 'MarkerSize',10)
ylabel('Velocity norm (m/s)')
xlabel('Time (hours)')
leg = legend('Cannonball SRP','Panel SRP');
set(leg, 'interpreter','latex')
set(gca,'FontSize',14)
grid on

%%
% load spire83_propagation_2018_11_7
figure(2)
% subplot(2,1,1)
plot(sod_uni/3600, theta, 'LineWidth',1)
hold on
plot(sod_sec_sp3log/3600, theta_sp3log*180/pi ,'k.', 'MarkerSize',10)
ylabel('$\theta$', 'interpreter','latex')
xlabel('Time (hours)')
title('Angle of velocity vector in body-frame')
leg = legend('Interpolated','Data');
set(leg, 'interpreter','latex')
set(gca,'FontSize',14)
grid on
%% EKF PLOTS
clc
clearvars
load est_msisFourier_srpcball
ys_res_postfit1 = ys_res_postfit;
rho_est1 = rho_est;
load est_msisFourier_srpcball_iter2
figure(1)
subplot(2,1,1)
plot(sod_sec_sp3/3600, vecnorm(ys_res_postfit1(1:3,:),2,1),'.', 'MarkerSize',10)
hold on
plot(sod_sec_sp3/3600, vecnorm(ys_res_postfit(1:3,:),2,1),'.', 'MarkerSize',10)
leg = legend('Iteration 1','Iteration 2');
set(leg, 'interpreter','latex')
ylabel('Position norm (m)')
title('Measurement residuals')
set(gca,'FontSize',16)
grid on

subplot(2,1,2)
plot(sod_sec_sp3/3600, vecnorm(ys_res_postfit1(4:6,:),2,1),'.', 'MarkerSize',10)
hold on
plot(sod_sec_sp3/3600, vecnorm(ys_res_postfit(4:6,:),2,1),'.', 'MarkerSize',10)
leg = legend('Iteration 1','Iteration 2');
set(leg, 'interpreter','latex')
ylabel('Velocity norm (m/s)')
xlabel('Time (hours)')
set(gca,'FontSize',16)
grid on


%%
figure(1)
load('spire83_ekf_jb08bias_iter2','rho_est','rho_nom','time_prop_utc_ekf')
load('spire83_ekf_msisbias','rho_nom')
plot(time_prop_utc_ekf/3600, rho_nom,'r.', 'MarkerSize',5)
hold on
plot(time_prop_utc_ekf/3600, rho_est,'k.', 'MarkerSize',5)
% load('spire83_rhohasdm','rho_est','rho_nom')
load('spire83_rhohasdm','rho_nom')
plot(time_prop_utc_ekf/3600, rho_nom,'b.', 'MarkerSize',5)
% plot(time_prop_utc_ekf/3600, rho_est,'g.', 'MarkerSize',5)
ylabel('Density ($kg/m^3$)','interpreter','latex')
xlabel('Time (hours)')
title('Estimated density')
% leg = legend('Estimated (Iteration 1)','Estimated (Iteration 2)','NRLMSISE-00','HASDM');
leg = legend('NRLMSISE-00','Estimate (MSIS00)','HASDM', 'Estimate (HASDM)');
set(leg, 'interpreter','latex')
set(gca,'FontSize',16)
grid on
% xlim([1 25])

%% error computations
clc
clearvars
load spire83_rhohasdm
rho_hasdm = rho_nom;
load spire83_ekf_jb08bias
ys_res_postfit1 = ys_res_postfit;
rho_est1 = rho_est;
load spire83_ekf_jb08bias_iter2

rho_err1 = abs(rho_hasdm - rho_est1)./rho_hasdm*100;
rho_err2 = abs(rho_hasdm - rho_est)./rho_hasdm*100;
rho_err0 = abs(rho_hasdm - rho_nom)./rho_hasdm*100;
rhom0 = mean(rho_err0)
rhorms0 = rms(rho_err0)
rhom1 = mean(rho_err1)
rhorms1 = rms(rho_err1)
rhom2 = mean(rho_err2)
rhorms2 = rms(rho_err2)

pos1 = vecnorm(ys_res_postfit1(1:3,:),2,1);
posm1 = mean(pos1)
posrms1 = rms(pos1)

vel1 = vecnorm(ys_res_postfit1(4:6,:),2,1);
velm1 = mean(vel1)
velrms1 = rms(vel1)

pos = vecnorm(ys_res_postfit(1:3,:),2,1);
posm2 = mean(pos)
posrms2 = rms(pos)

vel = vecnorm(ys_res_postfit(4:6,:),2,1);
velm2 = mean(vel)
velrms2 = rms(vel)