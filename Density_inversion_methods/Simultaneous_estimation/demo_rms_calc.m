clc
clearvars
%% Inputs
year_data = '2022';
month_data = '02';
day_data = '05';
%% Extract datafiles 
dir_data = '/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods/data/demo';
date_curr = strcat(year_data,'_',month_data,'_',day_data);
addpath(dir_data)
dir_struct = dir(dir_data);
dir_name = {dir_struct.name};
%% Spire corrected

dir_name1 = dir_name(contains(dir_name, 'spireResults_JB08'));
postfit_spire_corr_jb08 = [];
for ii = 1:numel(dir_name1)
    load(dir_name1{ii})
    if ~isempty(rho_est_trunc)
    postfit_spire_corr_jb08 = [postfit_spire_corr_jb08,ys_res_postfit];
    end
end

dir_name2 = dir_name(contains(dir_name, 'spireResults_MSIS00'));
postfit_spire_corr_msis = [];
for ii = 1:numel(dir_name2)
    load(dir_name2{ii})
    if ~isempty(rho_est_trunc)
    postfit_spire_corr_msis = [postfit_spire_corr_msis,ys_res_postfit];
    end
end

dir_name3 = dir_name(contains(dir_name, 'spireResults_WAMIPE'));
postfit_spire_corr_wamipe = [];
for ii = 1:numel(dir_name3)
    load(dir_name3{ii})
    if ~isempty(rho_est_trunc)
    postfit_spire_corr_wamipe = [postfit_spire_corr_wamipe,ys_res_postfit];
    end
end
%% Spire uncorrected

dir_name4 = dir_name(contains(dir_name, 'spireResults1_JB08'));
postfit_spire_jb08 = [];
for ii = 1:numel(dir_name4)
    load(dir_name4{ii})
    if ~isempty(rho_est_trunc)
    postfit_spire_jb08 = [postfit_spire_jb08,ys_res_postfit];
    end
end

dir_name5 = dir_name(contains(dir_name, 'spireResults1_MSIS00'));
postfit_spire_msis = [];
for ii = 1:numel(dir_name5)
    load(dir_name5{ii})
    if ~isempty(rho_est_trunc)
    postfit_spire_msis = [postfit_spire_msis,ys_res_postfit];
    end
end

dir_name6 = dir_name(contains(dir_name, 'spireResults1_WAMIPE'));
postfit_spire_wamipe = [];
for ii = 1:numel(dir_name6)
    load(dir_name6{ii})
    if ~isempty(rho_est_trunc)
    postfit_spire_wamipe = [postfit_spire_wamipe,ys_res_postfit];
    end
end

%% Starlink corrected

dir_name7 = dir_name(contains(dir_name, 'starlinkResults_JB08'));
postfit_starlink_corr_jb08 = [];
for ii = 1:numel(dir_name7)
    load(dir_name7{ii})
    if ~isempty(rho_est_trunc)
    postfit_starlink_corr_jb08 = [postfit_starlink_corr_jb08,ys_res_postfit];
    end
end

dir_name8 = dir_name(contains(dir_name, 'starlinkResults_MSIS00'));
postfit_starlink_corr_msis = [];
for ii = 1:numel(dir_name8)
    load(dir_name8{ii})
    if ~isempty(rho_est_trunc)
    postfit_starlink_corr_msis = [postfit_starlink_corr_msis,ys_res_postfit];
    end
end

dir_name9 = dir_name(contains(dir_name, 'starlinkResults_WAMIPE'));
postfit_starlink_corr_wamipe = [];
for ii = 1:numel(dir_name9)
    load(dir_name9{ii})
    if ~isempty(rho_est_trunc)
    postfit_starlink_corr_wamipe = [postfit_starlink_corr_wamipe,ys_res_postfit];
    end
end

%% Starlink uncorrected

dir_name10 = dir_name(contains(dir_name, 'starlinkResults1_JB08'));
postfit_starlink_jb08 = [];
for ii = 1:numel(dir_name10)
    load(dir_name10{ii})
    if ~isempty(rho_est_trunc)
    postfit_starlink_jb08 = [postfit_starlink_jb08,ys_res_postfit];
    end
end

dir_name11 = dir_name(contains(dir_name, 'starlinkResults1_MSIS00'));
postfit_starlink_msis = [];
for ii = 1:numel(dir_name11)
    load(dir_name11{ii})
    if ~isempty(rho_est_trunc)
    postfit_starlink_msis = [postfit_starlink_msis,ys_res_postfit];
    end
end
dir_name12 = dir_name(contains(dir_name, 'starlinkResults1_WAMIPE'));
postfit_starlink_wamipe = [];
for ii = 1:numel(dir_name12)
    load(dir_name12{ii})
    if ~isempty(rho_est_trunc)
    postfit_starlink_wamipe = [postfit_starlink_wamipe,ys_res_postfit];
    end
end

%% Plots
rms_values(1,1) = rms(vecnorm(postfit_spire_corr_jb08(1:3,:),2,1));
rms_values(1,2) = rms(vecnorm(postfit_starlink_corr_jb08(1:3,:),2,1));
rms_values(2,1) = rms(vecnorm(postfit_spire_jb08(1:3,:),2,1));
rms_values(2,2) = rms(vecnorm(postfit_starlink_jb08(1:3,:),2,1));

rms_values(3,1) = rms(vecnorm(postfit_spire_corr_msis(1:3,:),2,1));
rms_values(3,2) = rms(vecnorm(postfit_starlink_corr_msis(1:3,:),2,1));
rms_values(4,1) = rms(vecnorm(postfit_spire_msis(1:3,:),2,1));
rms_values(4,2) = rms(vecnorm(postfit_starlink_msis(1:3,:),2,1));

rms_values(5,1) = rms(vecnorm(postfit_spire_corr_wamipe(1:3,:),2,1));
rms_values(5,2) = rms(vecnorm(postfit_starlink_corr_wamipe(1:3,:),2,1));
rms_values(6,1) = rms(vecnorm(postfit_spire_wamipe(1:3,:),2,1));
rms_values(6,2) = rms(vecnorm(postfit_starlink_wamipe(1:3,:),2,1));

%%
X = categorical([{'Spire'}, {'Starlink'}]);
b = bar(X, rms_values);
ylabel('RMS (m)')
title('RMS values of postfit residuals')
legend('Corrected JB2008','JB2008','Corrected MSIS00','MSIS00','Corrected WAMIPE','WAMIPE')
set(gca, 'FontSize', 16)
grid on
set(gca,'YScale','log')
b(1).FaceColor = [1 0 0];
b(1).LineStyle = "-";
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(2).LineStyle = ":";
b(3).FaceColor = [0 1 1];
b(4).FaceColor = [0.3010 0.7450 0.9330];
b(4).LineStyle = ":";
b(5).FaceColor = [0 1 0];
b(6).FaceColor = [0.4660 0.6740 0.1880];
b(6).LineStyle = ":";