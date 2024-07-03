%% Monte Carlo cases
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/DataTools/mice/mice/src/mice/')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/DataTools/mice/mice/lib/' )
addpath('/Users/vishalray/GoogleDrive/Vishal/JB2K8/')
% addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Orbit - Processing real data')
clc
clear all

% addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Observability study/Results/OFF/Exponential')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Observability study/Results/BFF')



%% Case 4: {A00,A01,A40,A41}
% clear all
% noise_file = 'MC_noise_bff_mod';
% load(noise_file,'N_cases')
% str_est = [{'A00'}];%[{'A00'},{'A01'},{'A40'},{'A41'}];
% cd_file = 'Cd_gpm_dria_150F10';
% truth_file = 'Truth_gpm_dria_iner_bff';
% for ind_case = 1
%     ind_case
%     [X_state_err,Cd_err_rms,rms_res_f,rms_res_norm,del_X,Af_err_init,Af_err_est,P_up,Xf_true_all,Xf_name_all] =...
%         batch_cases(str_est,noise_file,cd_file,truth_file,ind_case,'mod');
%     mc_Xerr_est(:,ind_case) = X_state_err;
%     mc_poserr(:,ind_case) = norm(X_state_err(1:3));
%     mc_velerr(:,ind_case) = norm(X_state_err(4:6));
%     mc_cderr(:,ind_case) = Cd_err_rms;
%     mc_rms_res(:,ind_case) = rms_res_f;
%     mc_rms_res_norm(:,ind_case) = rms_res_norm;
%     mc_Xerr_init(:,ind_case) = del_X;
%     mc_Aferr_init(:,ind_case) = Af_err_init';
%     mc_Aferr_est(:,ind_case) = Af_err_est';
% end
% save('MC_bodf_case1_gpm_mod','P_up','Xf_true_all','Xf_name_all','-append')
% 
% clear all
% noise_file = 'MC_noise_bff_gpm';
% load(noise_file,'N_cases')
% str_est = [{'A00'},{'A20'}];%[{'A00'},{'A01'},{'A40'},{'A41'}];
% cd_file = 'Cd_gpm_dria_150F10';
% truth_file = 'Truth_gpm_dria_iner_bff';
% for ind_case = 1
%     ind_case
%     [X_state_err,Cd_err_rms,rms_res_f,rms_res_norm,del_X,Af_err_init,Af_err_est,P_up,Xf_true_all,Xf_name_all] =...
%         batch_cases(str_est,noise_file,cd_file,truth_file,ind_case,'mod');
%     mc_Xerr_est(:,ind_case) = X_state_err;
%     mc_poserr(:,ind_case) = norm(X_state_err(1:3));
%     mc_velerr(:,ind_case) = norm(X_state_err(4:6));
%     mc_cderr(:,ind_case) = Cd_err_rms;
%     mc_rms_res(:,ind_case) = rms_res_f;
%     mc_rms_res_norm(:,ind_case) = rms_res_norm;
%     mc_Xerr_init(:,ind_case) = del_X;
%     mc_Aferr_init(:,ind_case) = Af_err_init';
%     mc_Aferr_est(:,ind_case) = Af_err_est';
% end
% save('MC_bodf_case2_gpm_mod','P_up','Xf_true_all','Xf_name_all','-append')
% 
% clear all
% noise_file = 'MC_noise_bff_gpm';
% load(noise_file,'N_cases')
% str_est = [{'A00'},{'A20'},{'B10'}];%[{'A00'},{'A01'},{'A40'},{'A41'}];
% cd_file = 'Cd_gpm_dria_150F10';
% truth_file = 'Truth_gpm_dria_iner_bff';
% for ind_case = 1
%     ind_case
%     [X_state_err,Cd_err_rms,rms_res_f,rms_res_norm,del_X,Af_err_init,Af_err_est,P_up,Xf_true_all,Xf_name_all] =...
%         batch_cases(str_est,noise_file,cd_file,truth_file,ind_case,'mod');
%     mc_Xerr_est(:,ind_case) = X_state_err;
%     mc_poserr(:,ind_case) = norm(X_state_err(1:3));
%     mc_velerr(:,ind_case) = norm(X_state_err(4:6));
%     mc_cderr(:,ind_case) = Cd_err_rms;
%     mc_rms_res(:,ind_case) = rms_res_f;
%     mc_rms_res_norm(:,ind_case) = rms_res_norm;
%     mc_Xerr_init(:,ind_case) = del_X;
%     mc_Aferr_init(:,ind_case) = Af_err_init';
%     mc_Aferr_est(:,ind_case) = Af_err_est';
% end
% save('MC_bodf_case3_gpm_mod','P_up','Xf_true_all','Xf_name_all','-append')
% 
% clear all
% noise_file = 'MC_noise_bff_gpm';
% load(noise_file,'N_cases')
% str_est = [{'A00'},{'A20'},{'B10'},{'A40'}];%[{'A00'},{'A01'},{'A40'},{'A41'}];
% cd_file = 'Cd_gpm_dria_150F10';
% truth_file = 'Truth_gpm_dria_iner_bff';
% for ind_case = 1
%     ind_case
%     [X_state_err,Cd_err_rms,rms_res_f,rms_res_norm,del_X,Af_err_init,Af_err_est,P_up,Xf_true_all,Xf_name_all] =...
%         batch_cases(str_est,noise_file,cd_file,truth_file,ind_case,'mod');
%     mc_Xerr_est(:,ind_case) = X_state_err;
%     mc_poserr(:,ind_case) = norm(X_state_err(1:3));
%     mc_velerr(:,ind_case) = norm(X_state_err(4:6));
%     mc_cderr(:,ind_case) = Cd_err_rms;
%     mc_rms_res(:,ind_case) = rms_res_f;
%     mc_rms_res_norm(:,ind_case) = rms_res_norm;
%     mc_Xerr_init(:,ind_case) = del_X;
%     mc_Aferr_init(:,ind_case) = Af_err_init';
%     mc_Aferr_est(:,ind_case) = Af_err_est';
% end
% save('MC_bodf_case4_gpm_mod','P_up','Xf_true_all','Xf_name_all','-append')
% 
% 
% clear all
% noise_file = 'MC_noise_bff_gpm';
% load(noise_file,'N_cases')
% str_est = [{'A00'},{'A20'},{'B10'},{'A40'},{'A01'}];%[{'A00'},{'A01'},{'A40'},{'A41'}];
% cd_file = 'Cd_gpm_dria_150F10';
% truth_file = 'Truth_gpm_dria_iner_bff';
% for ind_case = 1
%     ind_case
%     [X_state_err,Cd_err_rms,rms_res_f,rms_res_norm,del_X,Af_err_init,Af_err_est,P_up,Xf_true_all,Xf_name_all] =...
%         batch_cases(str_est,noise_file,cd_file,truth_file,ind_case,'mod');
%     mc_Xerr_est(:,ind_case) = X_state_err;
%     mc_poserr(:,ind_case) = norm(X_state_err(1:3));
%     mc_velerr(:,ind_case) = norm(X_state_err(4:6));
%     mc_cderr(:,ind_case) = Cd_err_rms;
%     mc_rms_res(:,ind_case) = rms_res_f;
%     mc_rms_res_norm(:,ind_case) = rms_res_norm;
%     mc_Xerr_init(:,ind_case) = del_X;
%     mc_Aferr_init(:,ind_case) = Af_err_init';
%     mc_Aferr_est(:,ind_case) = Af_err_est';
% end
% save('MC_bodf_case5_gpm_mod','P_up','Xf_true_all','Xf_name_all','-append')
% 
% clear all
% noise_file = 'MC_noise_bff_gpm';
% load(noise_file,'N_cases')
% str_est = [{'A00'},{'A20'},{'B10'},{'A40'},{'A01'},{'A60'}];%[{'A00'},{'A01'},{'A40'},{'A41'}];
% cd_file = 'Cd_gpm_dria_150F10';
% truth_file = 'Truth_gpm_dria_iner_bff';
% for ind_case = 1
%     ind_case
%     [X_state_err,Cd_err_rms,rms_res_f,rms_res_norm,del_X,Af_err_init,Af_err_est,P_up,Xf_true_all,Xf_name_all] =...
%         batch_cases(str_est,noise_file,cd_file,truth_file,ind_case,'mod');
%     mc_Xerr_est(:,ind_case) = X_state_err;
%     mc_poserr(:,ind_case) = norm(X_state_err(1:3));
%     mc_velerr(:,ind_case) = norm(X_state_err(4:6));
%     mc_cderr(:,ind_case) = Cd_err_rms;
%     mc_rms_res(:,ind_case) = rms_res_f;
%     mc_rms_res_norm(:,ind_case) = rms_res_norm;
%     mc_Xerr_init(:,ind_case) = del_X;
%     mc_Aferr_init(:,ind_case) = Af_err_init';
%     mc_Aferr_est(:,ind_case) = Af_err_est';
% end
% save('MC_bodf_case6_gpm_mod','P_up','Xf_true_all','Xf_name_all','-append')
% 
% clear all
% noise_file = 'MC_noise_bff_gpm';
% load(noise_file,'N_cases')
% str_est = [{'A00'},{'A20'},{'B10'},{'A40'},{'A01'},{'A60'},{'A21'}];%[{'A00'},{'A01'},{'A40'},{'A41'}];
% cd_file = 'Cd_gpm_dria_150F10';
% truth_file = 'Truth_gpm_dria_iner_bff';
% for ind_case = 1
%     ind_case
%     [X_state_err,Cd_err_rms,rms_res_f,rms_res_norm,del_X,Af_err_init,Af_err_est,P_up,Xf_true_all,Xf_name_all] =...
%         batch_cases(str_est,noise_file,cd_file,truth_file,ind_case,'mod');
%     mc_Xerr_est(:,ind_case) = X_state_err;
%     mc_poserr(:,ind_case) = norm(X_state_err(1:3));
%     mc_velerr(:,ind_case) = norm(X_state_err(4:6));
%     mc_cderr(:,ind_case) = Cd_err_rms;
%     mc_rms_res(:,ind_case) = rms_res_f;
%     mc_rms_res_norm(:,ind_case) = rms_res_norm;
%     mc_Xerr_init(:,ind_case) = del_X;
%     mc_Aferr_init(:,ind_case) = Af_err_init';
%     mc_Aferr_est(:,ind_case) = Af_err_est';
% end
% save('MC_bodf_case7_gpm_mod','P_up','Xf_true_all','Xf_name_all','-append')
% 
% 
% clear all
% noise_file = 'MC_noise_bff_gpm';
% load(noise_file,'N_cases')
% str_est = [{'A00'},{'A20'},{'B10'},{'A40'},{'A01'},{'A60'},{'A21'},{'B11'}];%[{'A00'},{'A01'},{'A40'},{'A41'}];
% cd_file = 'Cd_gpm_dria_150F10';
% truth_file = 'Truth_gpm_dria_iner_bff';
% for ind_case = 1
%     ind_case
%     [X_state_err,Cd_err_rms,rms_res_f,rms_res_norm,del_X,Af_err_init,Af_err_est,P_up,Xf_true_all,Xf_name_all] =...
%         batch_cases(str_est,noise_file,cd_file,truth_file,ind_case,'mod');
%     mc_Xerr_est(:,ind_case) = X_state_err;
%     mc_poserr(:,ind_case) = norm(X_state_err(1:3));
%     mc_velerr(:,ind_case) = norm(X_state_err(4:6));
%     mc_cderr(:,ind_case) = Cd_err_rms;
%     mc_rms_res(:,ind_case) = rms_res_f;
%     mc_rms_res_norm(:,ind_case) = rms_res_norm;
%     mc_Xerr_init(:,ind_case) = del_X;
%     mc_Aferr_init(:,ind_case) = Af_err_init';
%     mc_Aferr_est(:,ind_case) = Af_err_est';
% end
% save('MC_bodf_case8_gpm_mod','P_up','Xf_true_all','Xf_name_all','-append')
% 
% clear all
% noise_file = 'MC_noise_bff_gpm';
% load(noise_file,'N_cases')
% str_est = [{'A00'},{'A20'},{'B10'},{'A40'},{'A01'},{'A60'},{'A21'},{'B11'},{'A23'}];%[{'A00'},{'A01'},{'A40'},{'A41'}];
% cd_file = 'Cd_gpm_dria_150F10';
% truth_file = 'Truth_gpm_dria_iner_bff';
% for ind_case = 1
%     ind_case
%     [X_state_err,Cd_err_rms,rms_res_f,rms_res_norm,del_X,Af_err_init,Af_err_est,P_up,Xf_true_all,Xf_name_all] =...
%         batch_cases(str_est,noise_file,cd_file,truth_file,ind_case,'mod');
%     mc_Xerr_est(:,ind_case) = X_state_err;
%     mc_poserr(:,ind_case) = norm(X_state_err(1:3));
%     mc_velerr(:,ind_case) = norm(X_state_err(4:6));
%     mc_cderr(:,ind_case) = Cd_err_rms;
%     mc_rms_res(:,ind_case) = rms_res_f;
%     mc_rms_res_norm(:,ind_case) = rms_res_norm;
%     mc_Xerr_init(:,ind_case) = del_X;
%     mc_Aferr_init(:,ind_case) = Af_err_init';
%     mc_Aferr_est(:,ind_case) = Af_err_est';
% end
% save('MC_bodf_case9_gpm_mod','P_up','Xf_true_all','Xf_name_all','-append')
% 
% % %% Case 5: {A00,A01,A40,A41,A02}
% % disp('Case 5')
% % clear all
% % load('MC_noise_bff_mod','N_cases')
% % str_est = [{'A00'},{'A01'},{'A40'},{'A41'},{'A02'}];
% % noise_file = 'MC_noise_bff_mod';
% % cd_file = 'Cd_cube_dria_150F10';
% % truth_file = 'Truth_cube_dria_iner_bff';
% % parfor ind_case = 1:N_cases
% %     ind_case
% %     [X_state_err,Cd_err_rms,rms_res_f,rms_res_norm,del_X,Af_err_init,Af_err_est] = batch_cases(str_est,noise_file,cd_file,truth_file,ind_case,'mod');
% %     mc_Xerr_est(:,ind_case) = X_state_err;
% %     mc_poserr(:,ind_case) = norm(X_state_err(1:3));
% %     mc_velerr(:,ind_case) = norm(X_state_err(4:6));
% %     mc_cderr(:,ind_case) = Cd_err_rms;
% %     mc_rms_res(:,ind_case) = rms_res_f;
% %     mc_rms_res_norm(:,ind_case) = rms_res_norm;
% %     mc_Xerr_init(:,ind_case) = del_X;
% %     mc_Aferr_init(:,ind_case) = Af_err_init';
% %     mc_Aferr_est(:,ind_case) = Af_err_est';
% % end
% % save('MC_bodf_case5_cube_mod')
% % 
% % %% Case 6: {A00,A01,A40,A02,A03}
% % disp('Case 6')
% % clear all
% % load('MC_noise_bff_mod','N_cases')
% % str_est = [{'A00'},{'A01'},{'A40'},{'A41'},{'A02'},{'A03'}];
% % noise_file = 'MC_noise_bff_mod';
% % cd_file = 'Cd_cube_dria_150F10';
% % truth_file = 'Truth_cube_dria_iner_bff';
% % parfor ind_case = 1:N_cases
% %     ind_case
% %     [X_state_err,Cd_err_rms,rms_res_f,rms_res_norm,del_X,Af_err_init,Af_err_est] = batch_cases(str_est,noise_file,cd_file,truth_file,ind_case,'mod');
% %     mc_Xerr_est(:,ind_case) = X_state_err;
% %     mc_poserr(:,ind_case) = norm(X_state_err(1:3));
% %     mc_velerr(:,ind_case) = norm(X_state_err(4:6));
% %     mc_cderr(:,ind_case) = Cd_err_rms;
% %     mc_rms_res(:,ind_case) = rms_res_f;
% %     mc_rms_res_norm(:,ind_case) = rms_res_norm;
% %     mc_Xerr_init(:,ind_case) = del_X;
% %     mc_Aferr_init(:,ind_case) = Af_err_init';
% %     mc_Aferr_est(:,ind_case) = Af_err_est';
% % end
% % save('MC_bodf_case6_cube_mod')
% % 
% 
%% Case 4: {A00,A01,A40,A41}
disp('Case 2')
clear all
load('MC_noise_bff_mod','N_cases')
str_est = [{'A00'}];
noise_file = 'MC_noise_bff_mod';
cd_file = 'Cd_cube_dria_150F10';
truth_file = 'Truth_cube_dria_iner_bff';
parfor ind_case = 1:N_cases
    ind_case
    [X_state_err,Cd_err_rms,rms_res_f,rms_res_norm,del_X,Af_err_init,Af_err_est] = batch_cases(str_est,noise_file,cd_file,truth_file,ind_case,'model');
    mc_Xerr_est(:,ind_case) = X_state_err;
    mc_poserr(:,ind_case) = norm(X_state_err(1:3));
    mc_velerr(:,ind_case) = norm(X_state_err(4:6));
    mc_cderr(:,ind_case) = Cd_err_rms;
    mc_rms_res(:,ind_case) = rms_res_f;
    mc_rms_res_norm(:,ind_case) = rms_res_norm;
    mc_Xerr_init(:,ind_case) = del_X;
    mc_Aferr_init(:,ind_case) = Af_err_init';
    mc_Aferr_est(:,ind_case) = Af_err_est';
end
% save('MC_bodf_case4_cube')


%% Function to run the filter
function [Xstate_err,Cd_err_rms,rms_res_f,rms_res_norm,del_X,Af_err_init,Af_err_est,P_up,Xf_true_all,Xf_name_all] =...
    batch_cases(str_est,noise_file,cd_file,truth_file,ind_case,str_case)
% load('Cd_sphere_dria_150F10','Af_total_mat','Af_mean','Bf_mean','Bf_std','Kl_mat')
load(cd_file,'Af_total_mat','Bf_total_mat','Cf_total_mat','Df_total_mat','Af_std','Bf_std','Cf_std','Df_std','Kl_mat','Af_Kl')
% load('Cd_gpm_dria_150F10','Af_total_mat','Af_mean','Af_std','Kl_mat','Af_Kl','Bf_Kl','Bf_total_mat','Bf_mean','Bf_std')
load(noise_file)
truth_model = truth_file;% 'Truth_sphere_dria';%'Truth_cube_dria_iner_bff';  %;
%% Flags
case_run = 'consider';                % 'truth', 'estimation', 'consider'
% N_size = size(Af_total_mat);
K_ind = 1; %N_size(end);               % true drag-coeff selection

flag_rho = 'Exp';                  % 'MSIS00', 'JB08', 'Exp'
flag_srp = 'Cball';                % Cball - cannonball, Three = three constants
estimated_coeff = [{'Cd'}];        % 'Cr','Cd'
flag_drag = 'Bod_orb';                 % body model: Bod, orbit model: Orb; Body orbit model: Bod_orb
bo_mod = 'BODF';                   % BODF : don't forget to change both otheta and ophi
Order_b = 30;
Order_o = 10;
del_X = del_X_mat(:,ind_case);
Pos_noise = Pos_noise_mat(:,:,ind_case);
Vel_noise = Vel_noise_mat(:,:,ind_case);
if strcmp(flag_drag,'Orb')
    Af_noise = Af_noise_mat(ind_case,:);
elseif strcmp(flag_drag,'Bod')
    Af_noise = Af_noise_mat(:,ind_case);
elseif strcmp(flag_drag,'Bod_orb')
    Af_noise = Af_noise_mat(:,:,ind_case);
    Bf_noise = Bf_noise_mat(:,:,ind_case);
    Cf_noise = Cf_noise_mat(:,:,ind_case);
    Df_noise = Df_noise_mat(:,:,ind_case);
end
run Main_obs_mc
% if strcmp(flag_drag,'Bod')
%     Af_err_init = Af_err_init';
%     Af_err_est = Af_err_est';
% elseif strcmp(flag_drag,'Bod_orb')
%     Af_err_init = 0;%Af_err_init';
%     Af_err_est = 0; %Af_err_est';
% end
rms_res_f = rms_res_post(:,end);
rms_res_norm = norm(rms_res_f./sqrt(diag(R_aug)))/sqrt(6);
end