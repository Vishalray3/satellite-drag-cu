%%
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Density inversion methods/Simultaneous estimation/error analysis')
clc
clearvars
Hmat = [300:100:800]; 
gravmat = [30:10:120];
arclen_mat = [10, 60, 60*20, 60*40, 60*90, 60*100];
hasdm_model1_mat = [{'2007_HASDM_300-375KM'}, {'2007_HASDM_400-475KM'}, {'2007_HASDM_500-575KM'}, {'2007_HASDM_600-675KM'}, {'2007_HASDM_700-825KM'}, {'2007_HASDM_700-825KM'}];
hasdm_model2_mat = [{'2007_HASDM_175-275KM'}, {'2007_HASDM_300-375KM'}, {'2007_HASDM_400-475KM'}, {'2007_HASDM_500-575KM'}, {'2007_HASDM_600-675KM'}, {'2007_HASDM_700-825KM'}];
parfor ii_ind = 1:numel(Hmat)
    ii_ind
    Halt_ind = Hmat(ii_ind);
    Hp_ind = Hmat(ii_ind) + 5;
    Ha_ind = Hmat(ii_ind) + 10;
    arclen = arclen_mat(ii_ind);
    truth_model = strcat('Truth_', num2str(Halt_ind), 'km_grav200');
    hasdm_model1 = hasdm_model1_mat{ii_ind};
    hasdm_model2 = hasdm_model2_mat{ii_ind};
    [rho_est_mean, rho_est_rms]= EDR_method_par(truth_model, hasdm_model1, hasdm_model2, Halt_ind, Hp_ind, Ha_ind, gravmat, arclen);
    rho_est_mean_array(ii_ind,:) = rho_est_mean;
    rho_est_rms_array(ii_ind,:) = rho_est_rms;
end
save('PODAccel_gravityresults_200')

function [rho_est_mean_temp, rho_est_rms_temp]= EDR_method_par(truth_model, hasdm_model1, hasdm_model2, Halt_ind, Hp_ind, Ha_ind, gravmat,arclen)
    for jj_ind = 1:numel(gravmat)
    gravdeg = gravmat(jj_ind);
    run Main_simulStudy
    rho_est_mean_temp(1,jj_ind) = rho_est_mean;
    rho_est_rms_temp(1,jj_ind) = rho_est_rms;
%     clearvars -except  rho_est_mean_temp  rho_est_rms_temp  Halt_ind  Hp_ind  Ha_ind  truth_model  gravmat jj_ind  hasdm_model1  hasdm_model2 arclen
    end
end
