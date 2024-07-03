clc
clearvars
load spire83_ekf_msisbias
load('spire83_propagation_2018_11_7', 'step','rho_hasdm','rho_msis')
time_step = step(1:10:end);

rho_hasdm_t = interp1(time_step, rho_hasdm, time_prop_utc_ekf);
rho_msis_t = interp1(time_step, rho_msis, time_prop_utc_ekf);


err_rho_est = rho_hasdm_t - rho_est;
mean_est = mean(err_rho_est);
rms_est = rms(err_rho_est);


err_rho_msis = rho_hasdm_t - rho_msis_t;
mean_msis = mean(err_rho_msis);
rms_msis = rms(err_rho_msis);
(mean_est-mean_msis)/mean_msis*100
(rms_est-rms_msis)/rms_msis*100
clearvars
load spire83_ekf_msisbias_iter2
load('spire83_propagation_2018_11_7', 'step','rho_hasdm','rho_msis')
time_step = step(1:10:end);

rho_hasdm_t = interp1(time_step, rho_hasdm, time_prop_utc_ekf);
rho_msis_t = interp1(time_step, rho_msis, time_prop_utc_ekf);


err_rho_est = rho_hasdm_t - rho_est;
mean_est = mean(err_rho_est);
rms_est = rms(err_rho_est);

err_rho_msis = rho_hasdm_t - rho_msis_t;
mean_msis = mean(err_rho_msis);
rms_msis = rms(err_rho_msis);
(mean_est-mean_msis)/mean_msis*100
(rms_est-rms_msis)/rms_msis*100


% plot(time_prop_utc_ekf, rho_hasdm_t,'.')
% hold on
% plot(time_prop_utc_ekf, rho_msis_t,'.')
% plot(time_prop_utc_ekf, rho_est,'.')
% 



