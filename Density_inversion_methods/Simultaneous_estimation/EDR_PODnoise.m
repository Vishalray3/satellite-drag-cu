%%
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Density inversion methods/Simultaneous estimation/error analysis/EDR')
clc
clearvars
rng('default')
Hmat = [300:100:800];

arclen_mat = [10, 60, 10*60, [30:30:1200]*60];
roll_mat = [0];%:1:10];
pitch_mat = [0];%:1:10];
hasdm_model1_mat = [{'2007_HASDM_300-375KM'}, {'2007_HASDM_400-475KM'}, {'2007_HASDM_500-575KM'}, {'2007_HASDM_600-675KM'}, {'2007_HASDM_700-825KM'}, {'2007_HASDM_700-825KM'}];
hasdm_model2_mat = [{'2007_HASDM_175-275KM'}, {'2007_HASDM_300-375KM'}, {'2007_HASDM_400-475KM'}, {'2007_HASDM_500-575KM'}, {'2007_HASDM_600-675KM'}, {'2007_HASDM_700-825KM'}];
% hasdm_model1_mat = [ {'2007_HASDM_700-825KM'}];
% hasdm_model2_mat = [{'2007_HASDM_700-825KM'}];


% case_run = 'Truth';
% sigma_pos = 0.1;
% sigma_vel = 1e-4;
% parfor ii_ind = 1:numel(Hmat)
%     ii_ind
%     Halt_ind = Hmat(ii_ind);
%     Hp_ind = Hmat(ii_ind) + 5;
%     Ha_ind = Hmat(ii_ind) + 10;
%     truth_model = strcat('Truth', num2str(Halt_ind), 'km_solarmin');
%     hasdm_model1 = hasdm_model1_mat{ii_ind};
%     hasdm_model2 = hasdm_model2_mat{ii_ind};
%     run_truth(truth_model, hasdm_model1, hasdm_model2, Halt_ind, Hp_ind, Ha_ind,...
%     arclen_mat, pitch_mat, roll_mat, case_run);
% end



case_run = 'Estimation';
sigma_pos = 0.1;
sigma_vel = 1e-4;
parfor ii_ind = 1:numel(Hmat)
    ii_ind
    Halt_ind = Hmat(ii_ind);
    Hp_ind = Hmat(ii_ind) + 5;
    Ha_ind = Hmat(ii_ind) + 10;
    truth_model = strcat('Truth', num2str(Halt_ind), 'km_solarmin');
    hasdm_model1 = hasdm_model1_mat{ii_ind};
    hasdm_model2 = hasdm_model2_mat{ii_ind};
    [rho_est_mean, rho_est_rms, rho_est_mean_avg, rho_est_rms_avg]= EDR_method_par(truth_model, hasdm_model1, hasdm_model2, Halt_ind, Hp_ind, Ha_ind,...
        arclen_mat,pitch_mat, roll_mat, case_run, sigma_pos, sigma_vel);
    rho_est_mean_array_10(ii_ind,:) = rho_est_mean;
    rho_est_rms_array_10(ii_ind,:) = rho_est_rms;
    rho_est_mean_avg_array_10(ii_ind,:) = rho_est_mean_avg;
    rho_est_rms_avg_array_10(ii_ind,:) = rho_est_rms_avg;
    rho_est_rms
%     rho_est_mat(ii_ind, :) = rho_est;
%     rho_nom_mat = rho_nom;
end

case_run = 'Estimation';
sigma_pos = 0.01;
sigma_vel = 1e-5;
parfor ii_ind = 1:numel(Hmat)
    ii_ind
    Halt_ind = Hmat(ii_ind);
    Hp_ind = Hmat(ii_ind) + 5;
    Ha_ind = Hmat(ii_ind) + 10;
    truth_model = strcat('Truth', num2str(Halt_ind), 'km_solarmin');
    hasdm_model1 = hasdm_model1_mat{ii_ind};
    hasdm_model2 = hasdm_model2_mat{ii_ind};
    [rho_est_mean, rho_est_rms, rho_est_mean_avg, rho_est_rms_avg]= EDR_method_par(truth_model, hasdm_model1, hasdm_model2, Halt_ind, Hp_ind, Ha_ind,...
        arclen_mat,pitch_mat, roll_mat, case_run, sigma_pos, sigma_vel);
    rho_est_mean_array_1(ii_ind,:) = rho_est_mean;
    rho_est_rms_array_1(ii_ind,:) = rho_est_rms;
    rho_est_mean_avg_array_1(ii_ind,:) = rho_est_mean_avg;
    rho_est_rms_avg_array_1(ii_ind,:) = rho_est_rms_avg;
%     rho_est_mat(ii_ind, :) = rho_est;
%     rho_nom_mat = rho_nom;
end

save('PODnoise_results_solarmin')



function [rho_est_mean_temp, rho_est_rms_temp, rho_est_avg_mean, rho_est_avg_rms]= EDR_method_par(truth_model, hasdm_model1, hasdm_model2,...
    Halt_ind, Hp_ind, Ha_ind, arclen, pitch_mat, roll_mat, case_run, sigma_pos_new, sigma_vel_new)
        run Main_simulStudy
        rho_est_mean_temp = rho_est_mean;
        rho_est_rms_temp = rho_est_rms;
        rho_est_avg_mean = rho_est_mean_avg;
        rho_est_avg_rms = rho_est_rms_avg;
end

function [ ]= run_truth(truth_model, hasdm_model1, hasdm_model2, Halt_ind, Hp_ind, Ha_ind,...
    arclen, pitch_mat, roll_mat, case_run)
        run Main_simulStudy
        save(truth_model)
end