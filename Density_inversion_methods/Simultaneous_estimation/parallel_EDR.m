%% parallel filters
clc
clearvars
% addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Density inversion methods/Simultaneous estimation/error analysis')
% addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Density inversion methods')
roll_mat = [0:1:10];
pitch_mat = [0:1:10];
Nr = numel(roll_mat);

for ii = 1:Nr
    ii
%     roll = normrnd(0,roll_mat(ii),[1,8641]);
    for jj = 1:Nr
%         pitch = normrnd(0,pitch_mat(jj),[1,8641]);
        [rho_est_mean, rho_est_rms] = parallel_filter(roll_mat(ii), pitch_mat(jj), pitch_mat(jj));
        rho_est_mat(ii,jj,:) = rho_est;
        rho_mean_mat(ii,jj)  = rho_est_mean;
        rho_rms_mat(ii,jj)   = rho_est_rms;
        Cd_est_mat(ii,jj,:)  = Cd_est;
        Cd_mean_mat(ii,jj,:) = Cd_est_mean;
        Cd_rms_mat(ii,jj,:)  = Cd_est_rms;
    end
end

save('EDR_300km_nadirerror')


function [rho_est_mean, rho_est_rms] = parallel_filter(roll, pitch, yaw)

run Main_simulStudy
parameters.roll = roll;
parameters.pitch = pitch;
parameters.yaw = yaw;
run EDR_method

end