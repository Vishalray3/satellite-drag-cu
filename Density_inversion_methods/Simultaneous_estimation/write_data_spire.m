%% Post-processing
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods')
clc
clearvars
%% Inputs
year_data = '2019';
month_data = '12';
day_data = '01';
wam_filename = 'WAM_den_20220205';
%% Extract datafiles 
date_curr = strcat(year_data,'_',month_data,'_',day_data);
dir_data = strcat('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods/data/spire/',date_curr);
addpath(dir_data)
dir_struct = dir(dir_data);
dir_name = {dir_struct.name};
dir_name = dir_name(contains(dir_name, 'spireResults_JB08'));
% dir_name = dir_name(contains(dir_name, 'starlinkResults_MSIS00'));
sat_ID_mat = extractBetween(dir_name, 'FM', '_2019');
% sat_ID_mat = extractBetween(dir_name, 'MSIS00_', '_2022');
sat_id_num = str2double(sat_ID_mat);
data_mat = [];

%% Read WAM-IPE data (if needed)
% F_wamipe = wamipe_interpolant(wam_filename, str2double(year_data), str2double(month_data), str2double(day_data));
% jd_ref = 86400*GREGORIANtoJD_vector(1950,1,0) + 3600*0 + 60*0 + 0;
%% Write data
for ii = 1:numel(dir_name)
    ii
    sat_ID = sat_id_num(ii); 
    load(dir_name{ii})
    sat_ID = sat_ID*ones(1,numel(rho_est_trunc));
    rho_wamipe = zeros(1, numel(rho_est_trunc));
    data_id = [sat_ID', rho_est_trunc', jdutc_trunc',lat_trunc', long_trunc', alt_trunc'/1e3, rho_nom_trunc', rho_hasdm_trunc', rho_wamipe'];
    data_mat = [data_mat;data_id];
end
data_tab = array2table(data_mat);
data_tab.Properties.VariableNames = {'Sat ID','Estimated Density (kg/m3)','UTC Julian date (sec)','Latitude (deg)','Longitude (deg)','Altitude (km)','JB2008', 'HASDM','WAM-IPE'};
writetable(data_tab, 'KayhanDensities_2019_12_01.csv')

