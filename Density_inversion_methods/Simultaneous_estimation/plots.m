
%% Plots

%%
figure()
plot(time_prop_utc/3600,rho_true-rho_est, 'k','LineWidth', 1)
hold on
plot(time_prop_utc/3600,Xs_est(7,:)+ Xs_est(9,:), '--r','LineWidth', 1)
xlabel('Time (hours)')
ylabel('Correction')
legend('Truth', 'Estimate')
title('Density correction')
set(gca,'FontSize',18)
grid on
%%
figure()
plot(time_prop_utc/3600,rho_true, 'k','LineWidth', 2)
hold on
plot(time_prop_utc/3600,rho_est, '--r','LineWidth', 1)
plot(time_prop_utc/3600,X_iter(3).xs_est+Xs_est(9,:)+rho_est, 'g','LineWidth', 1)
xlabel('Time (hours)')
ylabel('Density')
legend('Truth (MSIS)', 'Filter (JB08)', 'Estimate')
title('Density estimate')
set(gca,'FontSize',18)
grid on
xlim([0 23])
%%
load Case3_jb08
figure()
plot(time_prop_utc/3600,rho_true-rho_est, 'k','LineWidth', 2)
hold on
plot(time_prop_utc/3600,Xs_est(9,:)+Xs_est(7,:), '--r','LineWidth', 1)
load Case3_jb08_cderr
plot(time_prop_utc/3600,Xs_est(9,:)+Xs_est(7,:), '-.b','LineWidth', 1)
% plot(time_prop_utc/3600,X_iter(3).xs_est+Xs_est(9,:)+rho_est, 'g','LineWidth', 1)
xlabel('Time (hours)')
ylabel('Correction')
legend('Truth', 'Estimate: Cd known', 'Estimate: Cd biased')
title('Density correction')
set(gca,'FontSize',18)
grid on

%% 
load Case2_gpm_AllFourierEst_Cdest.mat
Cd_calc = Cd_est;

figure()
plot(time_prop_utc/3600,rho_true-rho_est, 'k','LineWidth', 2)
hold on
load Case2_gpm_AllFourierEst.mat
plot(time_prop_utc/3600,Xs_est(9,:)+Xs_est(7,:), 'r','LineWidth', 1)
load Case2_gpm_AllFourierEst_secondIter.mat
plot(time_prop_utc/3600,Xs_est(9,:)+Xs_est(7,:), 'g','LineWidth', 1)
% plot(time_prop_utc/3600,X_iter(3).xs_est+Xs_est(9,:)+rho_est, 'g','LineWidth', 1)
xlabel('Time (hours)')
ylabel('Density error')
legend('Truth', 'Cd bias', 'Bias estimated')
title('Estimated density error')
set(gca,'FontSize',18)
grid on
xlim([1 23])

%%

figure()
load Truth_msis_circ_gpm_tseries
plot(time_prop_utc/3600,Cd_true/2, 'k','LineWidth', 1)
hold on
load Truth_jb08_circ_gpm_tseries_kl
plot(time_prop_utc/3600,Cd_true/2, '--r','LineWidth', 1)
xlabel('Time (hours)')
ylabel('Drag-coefficient')
legend('Cd true', 'Cd filter')
title('Drag-coefficient variation')
set(gca,'FontSize',18)
grid on
xlim([1 10])
% ylim([0.3 2.5])

%%

figure()
load Truth_msis_circ_gpm_tseries
subplot(2,1,1)
plot(time_prop_utc/3600,Af_tseries(1,:), 'k','LineWidth', 1)
hold on
load Truth_jb08_circ_gpm_tseries_kl
plot(time_prop_utc/3600,Af_tseries(1,:), '--r','LineWidth', 1)
ylabel('Order 0')
legend('A0 true', 'A0 filter')
title('Fourier-coefficient variation')
set(gca,'FontSize',18)
grid on
xlim([1 10])

subplot(2,1,2)
load Truth_msis_circ_gpm_tseries
plot(time_prop_utc/3600,Af_tseries(3,:), 'k','LineWidth', 1)
hold on
load Truth_jb08_circ_gpm_tseries_kl
plot(time_prop_utc/3600,Af_tseries(3,:), '--r','LineWidth', 1)
xlabel('Time (hours)')
ylabel('Order 2')
legend('A2 true', 'A2 filter')
set(gca,'FontSize',18)
grid on
xlim([1 10])
% ylim([0.3 2.5])

%% 
load Truth_msis_circ_gpm_tseries_day302_2003
Cd_true1 = Cd_true;
figure()
load Truth_jb08_circ_gpm_tseries_day302_2003_kl
plot(time_prop_utc/3600,Cd_true1/2-Cd_true/2, 'k','LineWidth', 1)
hold on
load Case5_2003_ekf_iter1
plot(time_prop_utc/3600,Cd_true/2-Cd_est/2, '--r','LineWidth', 1)
load Case5_2003_ekf_iter2
plot(time_prop_utc/3600,Cd_true/2-Cd_est/2, '-.b','LineWidth', 1)
% load Case4_allFourier_kl_thirdIter
% plot(time_prop_utc/3600,Cd_true/2-Cd_est/2, 'g','LineWidth', 1)
xlabel('Time (hours)')
ylabel('Error')
legend('Iteration 0', 'Iteration 1', 'Iteration 2', 'Iteration 3')
title('Drag-coefficient error')
set(gca,'FontSize',18)
grid on
xlim([1 10])

%% 

figure()
load Case4_allFourier_kl
plot(time_prop_utc/3600,rho_true-rho_est, 'k','LineWidth', 1)
hold on
plot(time_prop_utc/3600,Xs_est(9,:)+Xs_est(7,:), '--r','LineWidth', 1)
load Case4_allFourier_kl_secondIter
plot(time_prop_utc/3600,Xs_est(9,:)+Xs_est(7,:), '-.g','LineWidth', 1)
load Case4_allFourier_kl_thirdIter
plot(time_prop_utc/3600,Xs_est(9,:)+Xs_est(7,:), 'm','LineWidth', 1)
xlabel('Time (hours)')
ylabel('Error')
legend('Truth', 'Iteration 1', 'Iteration 2', 'Iteration 3')
title('Density correction')
set(gca,'FontSize',18)
grid on
ylim([-3e-12 3.5e-12])

%% 
figure()
load('Case5_2003_ekf_iter1', 'time_prop_utc', 'rho_true','rho_est', 'Xs_est')
plot(time_prop_utc/3600,rho_true-rho_est, 'k','LineWidth', 1)
hold on
plot(time_prop_utc/3600,Xs_est(9,:)+Xs_est(7,:), 'c','LineWidth', 1)
load('Case5_2003_ekf_iter2',  'Xs_est')
plot(time_prop_utc/3600,Xs_est(9,:)+Xs_est(7,:), 'r','LineWidth', 1)
xlabel('Time (hours)')
ylabel('Error')
legend('Truth', 'Iteration 1', 'Iteration 2')
title('Density correction')
set(gca,'FontSize',18)
grid on
% ylim([-3e-12 3.5e-12])
%%
load('Case5_2003_ekf_iter1', 'time_prop_utc', 'rho_true','rho_est', 'Xs_est')
figure()
plot(time_prop_utc/3600,rho_true, 'k','LineWidth', 2)
hold on
plot(time_prop_utc/3600,rho_est, 'c','LineWidth', 2)
load('Case5_2003_ekf_iter2',  'Xs_est')
plot(time_prop_utc/3600,rho_est+Xs_est(9,:)+Xs_est(7,:), 'r','LineWidth', 1)
% plot(time_prop_utc/3600,X_iter(3).xs_est+Xs_est(9,:)+rho_est, 'g','LineWidth', 1)
xlabel('Time (hours)')
ylabel('Density')
legend('Truth (MSIS00)', 'Filter (JB08)', 'Estimate')
title('Estimated density')
set(gca,'FontSize',18)
grid on
%%
clearvars
clc
load Case4_allFourier_kl_thirdIter
Cd_err = Cd_true - Cd_est;
mean(Cd_err)
rms(Cd_err)
mean(error_dens(1:7000))
rms(error_dens(1:7000))