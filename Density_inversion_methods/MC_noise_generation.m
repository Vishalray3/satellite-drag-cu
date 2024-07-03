%% Generating noise
clc
clearvars
% load('Cd_sphere_dria_150F10','Af_std')
load('Cd_gpm_dria_150F10','Af_std','Bf_std','Cf_std','Df_std')
N_cases = 200; % Monte Carlo cases
del_r = [10;10;10];
del_v = [0.01;0.01;0.01];
del_X_mat = [normrnd(0,del_r(1),[3,N_cases]);normrnd(0,del_v(1),[3,N_cases])];


sig_meas(1) = 1.5;
sig_meas(2) = 5e-3;
sigma_pos = sig_meas(1);
sigma_vel = sig_meas(2);
Pos_noise_mat = normrnd(0,sigma_pos,[3,8641,N_cases]);
Vel_noise_mat = normrnd(0,sigma_vel,[3,8641,N_cases]);

%%
Af_noise_mat = normrnd(0,repmat(Af_std,1,1,N_cases)); % normrnd(0,repmat(Af_std' ,N_cases,1)); for OFF!?>
Bf_noise_mat = normrnd(0,repmat(Bf_std,1,1,N_cases)); 
Cf_noise_mat = normrnd(0,repmat(Cf_std,1,1,N_cases)); 
Df_noise_mat = normrnd(0,repmat(Df_std,1,1,N_cases)); 
save('MC_noise_bff_gpm')