% close all
clearvars
load('spireResults_JB08_FM099_2022_02_05_new')


%%
rho_true = rho_hasdm_trunc;
rho_est = rho_est_trunc;
rho_nom = rho_nom_trunc;

time_diff = diff(time_prop_trunc);
time_ind = find(time_diff > 100);
time_ind = [0, time_ind, numel(time_prop_trunc)];

ind_start = 1;
ind_end = 0;

for mm = 1:numel(time_ind)-1
    rho_true_eff(mm) = mean(rho_true(time_ind(mm)+ind_start:time_ind(mm+1)-ind_end));
    rho_eff(mm) = mean(rho_est(time_ind(mm)+ind_start:time_ind(mm+1)-ind_end));
    rho_nom_eff = mean(rho_nom(time_ind(mm)+ind_start:time_ind(mm+1)-ind_end));
    rho_est_error_avg(mm) = (rho_true_eff(mm)-rho_eff(mm))./rho_true_eff(mm)*100;
    rho_nom_error_avg(mm) = (rho_true_eff(mm)-rho_nom_eff)./rho_true_eff(mm)*100;
end

rho_est_mean_avg = mean(rho_est_error_avg);
rho_est_rms_avg  = rms(rho_est_error_avg);
% figure()
plot(abs(rho_est_error_avg))
hold on
plot(abs(rho_nom_error_avg))

%%
close all
plot(sod_trunc, rho_hasdm_trunc, '.')
hold on
plot(sod_trunc, rho_est_trunc, '.')
plot(sod_trunc, rho_nom_trunc, '.')


