%%
% addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Density inversion methods/Simultaneous estimation/error analysis/EDR')
clc
clearvars
rng('default')
Hmat = [800,700,800,700,700,800,600,600,700,800,600,700,800,600,500,800,500,800,600,500,700];

arclen_mat = [10, 60, 10*60, [30:30:1200]*60];

hasdm_model1_mat = [{'2007_HASDM_700-375KM'}, {'2007_HASDM_400-475KM'}, {'2007_HASDM_500-575KM'}, {'2007_HASDM_600-675KM'}, {'2007_HASDM_700-825KM'}, {'2007_HASDM_700-825KM'}; ...
    {'2003_HASDM_300-375KM'}, {'2003_HASDM_400-475KM'}, {'2003_HASDM_500-575KM'}, {'2003_HASDM_600-675KM'}, {'2003_HASDM_700-825KM'}, {'2003_HASDM_700-825KM'}];

hasdm_model2_mat = [{'2007_HASDM_175-275KM'}, {'2007_HASDM_300-375KM'}, {'2007_HASDM_400-475KM'}, {'2007_HASDM_500-575KM'}, {'2007_HASDM_600-675KM'}, {'2007_HASDM_700-825KM'}; ...
    {'2003_HASDM_175-275KM'}, {'2003_HASDM_300-375KM'}, {'2003_HASDM_400-475KM'}, {'2003_HASDM_500-575KM'}, {'2003_HASDM_600-675KM'}, {'2003_HASDM_700-825KM'}];
sigma_pos_mat = [0.01, 0.05, 0.1, 0.25, 0.5];

area_mat = [0.75,0.5,1,0.75,1,2.5,0.5,0.75,2.5,5,1,5,0.5,2.5,0.5,0.75,0.75,1,5,1,0.5];   % actual areas: 0.05, 0.5, 5, 50 - > amr: 0.01, 0.1, 1, 10

solar_level = [{'solarmin'},{'solarmin'},{'solarmin'},{'solarmin'},{'solarmin'},{'solarmin'},{'solarmin'},{'solarmin'},{'solarmin'},{'solarmin'},{'solarmin'},{'solarmin'},{'solarmax'},{'solarmin'},{'solarmin'},{'solarmax'},{'solarmin'},{'solarmax'},{'solarmin'},{'solarmin'},{'solarmax'},{'solarmin'}];

case_run = 'EDR';

N_h = numel(Hmat);
N_a = numel(area_mat);
N_s = numel(sigma_pos_mat);
N_sol = numel(solar_level);
for ii_ind = 1:N_h
    ii_ind
    Halt_ind = Hmat(ii_ind);
    Hp_ind = Hmat(ii_ind) + 5;
    Ha_ind = Hmat(ii_ind) + 10;
    if strcmp(solar_level{ii_ind}, 'solarmin')
        kk = 1;
    elseif strcmp(solar_level{ii_ind}, 'solarmax')
        kk = 2;
    end
    [rho_est_mean, rho_est_rms, rho_est_mean_avg, rho_est_rms_avg, Psig, arcmin, arc10, arc20]...
        = main_for(Halt_ind, Hp_ind, Ha_ind, arclen_mat, case_run, area_mat(ii_ind), sigma_pos_mat, solar_level{ii_ind}, hasdm_model1_mat(kk,ii_ind), hasdm_model2_mat(kk,ii_ind));
    
    rho_est_mean_array(ii_ind,:) = rho_est_mean;
    rho_est_rms_array(ii_ind,:) = rho_est_rms;
    rho_est_mean_avg_array(ii_ind,:) = rho_est_mean_avg;
    rho_est_rms_avg_array(ii_ind,:) = rho_est_rms_avg;
    Psig_array(ii_ind,:) = Psig;
    arcmin_array(ii_ind,:) = arcmin;
    arc10_array(ii_ind,:) = arc10;
    arc20_array(ii_ind,:) = arc20;
end


save('EDRnoise_snr2')


function [rho_est_mean_array, rho_est_rms_array, rho_est_mean_avg_array, rho_est_rms_avg_array, Psig_array, arcmin_array, arc10_array, arc20_array]...
    = main_for(Halt_ind, Hp_ind, Ha_ind, arclen_mat, case_run, area_mat, sigma_pos_mat, solar_level, hasdm_model1_mat, hasdm_model2_mat)

N_s = numel(sigma_pos_mat);

area = area_mat;

sol_lev = solar_level;

hasdm_model1 = hasdm_model1_mat;
hasdm_model2 = hasdm_model2_mat;
try
    [X_true_aug, parameters_truth, rho_true] = run_truth( hasdm_model1, hasdm_model2, ...
        Halt_ind, Hp_ind, Ha_ind, arclen_mat, area, sol_lev);
    
    for kk_ind = 1: N_s
        sigma_pos = sigma_pos_mat(kk_ind);
        sigma_vel = sigma_pos*1e-3;
        
        [rho_est_mean, rho_est_rms, rho_est_mean_avg, rho_est_rms_avg, P_signal, arc_min, arc_10, arc_20]= EDR_method_par(hasdm_model1, hasdm_model2, ...
            Halt_ind, Hp_ind, Ha_ind, arclen_mat, case_run, sigma_pos, sigma_vel, area, sol_lev, X_true_aug, parameters_truth, rho_true);
        rho_est_mean_array(kk_ind) = rho_est_mean;
        rho_est_rms_array(kk_ind) = rho_est_rms;
        rho_est_mean_avg_array(kk_ind,:) = rho_est_mean_avg;
        rho_est_rms_avg_array(kk_ind,:) = rho_est_rms_avg;
        Psig_array(kk_ind) = P_signal;
        arcmin_array(kk_ind) = arc_min;
        arc10_array(kk_ind) = arc_10;
        arc20_array(kk_ind) = arc_20;
    end
catch
    rho_est_mean_array(1,:) = NaN;
    rho_est_rms_array(1,:) = NaN;
    rho_est_mean_avg_array(1,:) = NaN;
    rho_est_rms_avg_array(1,:) = NaN;
    Psig_array(1,:) = NaN;
    arcmin_array(1,:) = NaN;
    arc10_array(1,:) = NaN;
    arc20_array(1,:) = NaN;
end
end



function [rho_est_mean_temp, rho_est_rms_temp, rho_est_avg_mean, rho_est_avg_rms, P_signal, arc_min, arc_10, arc_20]= EDR_method_par(hasdm_model1, hasdm_model2,...
    Halt_ind, Hp_ind, Ha_ind, arclen, case_run, sigma_pos_new, sigma_vel_new, area_new, sol_lev, X_true_aug, parameters_truth, rho_true)
run Main_simulStudy
rho_est_mean_temp = rho_est_mean;
rho_est_rms_temp = rho_est_rms;
rho_est_avg_mean = rho_mean_avg_min;
rho_est_avg_rms = rho_rms_avg_min;

end

function [X_true_aug, parameters_truth, rho_true]= run_truth( hasdm_model1, hasdm_model2, Halt_ind, Hp_ind, Ha_ind,...
    arclen, area_new, sol_lev)
case_run = 'Truth';
run Main_simulStudy
parameters_truth = parameters;
rho_true = rho;
end