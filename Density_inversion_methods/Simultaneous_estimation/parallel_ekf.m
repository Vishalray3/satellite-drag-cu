%% parallel filters
clc
clearvars


Hp_ind = 805;
Ha_ind = 810;
Halt_ind = 800;
arclen_mat = [60*100];
roll_mat = [0:1:10];
pitch_mat = [0:1:10];
Nr = numel(roll_mat);
hasdm_model1_mat = '2007_HASDM_700-825KM';
hasdm_model2_mat = '2007_HASDM_700-825KM';

parfor ii = 1:Nr
    ii
%     roll = normrnd(0,roll_mat(ii),[1,8641]);
    for jj = 1:Nr
%         pitch = normrnd(0,pitch_mat(jj),[1,8641]);
        [ rho_est_mean, rho_est_rms,  rho_est_mean_avg, rho_est_rms_avg] = parallel_filter(roll_mat(ii), pitch_mat(jj),...
            arclen_mat, hasdm_model1_mat, hasdm_model2_mat, 'Truth800km_erp', Ha_ind, Hp_ind, Halt_ind);
        rho_mean_mat(ii,jj)  = rho_est_mean;
        rho_rms_mat(ii,jj)   = rho_est_rms;
        rho_mean_avg_mat(ii,jj)  = rho_est_mean_avg;
        rho_rms_avg_mat(ii,jj)   = rho_est_rms_avg;
    end
end

save('EDR_800km_nadirerror')


function [rho_est_mean, rho_est_rms, rho_est_mean_avg, rho_est_rms_avg] = parallel_filter(roll_ind, pitch_ind, arclen, hasdm_model1, hasdm_model2, ...
    truth_model, Ha_ind, Hp_ind, Halt_ind)

run EDR_method

end