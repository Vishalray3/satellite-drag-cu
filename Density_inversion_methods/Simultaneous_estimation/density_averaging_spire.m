close all
clc
clearvars
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/main_working_folder_data/spire/2022_02_05')
load('spireResults_JB08_FM115_2022_02_05_new')
%%
%% Orbit period
r_circ = a_sma;
mu_e = 398600435436096;
v_circ = sqrt(mu_e/r_circ);
n_mean = v_circ/r_circ;
T_orb = 2*pi/n_mean;
%%
arclen = 60*60; 

rho_true = rho_hasdm_trunc;
rho_est = rho_est_trunc;% Cd_est(ind_vec).*rho_est_trunc/0.31;
rho_nom = rho_nom_trunc;

time_diff = diff(time_prop_trunc);
time_ind = find(time_diff > 100);
time_ind = [1, time_ind, numel(time_prop_trunc)];

ind_start = 1;
ind_end = 0;

time_step = time_prop_trunc(2) - time_prop_trunc(1);

[rho_true_results, rho_nom_results, rho_est_results, alt_results] = den_avg_spire(rho_true, rho_est, rho_nom, arclen, time_step, time_ind, alt_trunc);

% remove outliers
[~, ind_out] = rmoutliers(rho_est_results./rho_nom_results);

ind_retain = ~ind_out;


rho_err_est = (rho_est_results - rho_true_results)./rho_true_results*100;
rho_err_nom = (rho_nom_results - rho_true_results)./rho_true_results*100;

rho_est_mean = mean(rho_err_est(ind_retain))
rho_est_rms  = rms(rho_err_est(ind_retain))

rho_nom_mean = mean(rho_err_nom(ind_retain))
rho_nom_rms  = rms(rho_err_nom(ind_retain))

% figure()
% plot(abs(rho_err_est))
% hold on
% plot(abs(rho_err_nom))

plot(abs(rho_err_est(ind_retain)))
hold on
plot(abs(rho_err_nom(ind_retain)))

%%
% % close all
% plot(sod_trunc, rho_hasdm_trunc, '.')
% hold on
% plot(sod_trunc, rho_est_trunc, '.')
% plot(sod_trunc, rho_nom_trunc, '.')


