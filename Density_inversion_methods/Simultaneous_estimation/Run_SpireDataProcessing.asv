%% Live demo script
clc
clearvars
tic
restoredefaultpath
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods/data/HASDM_data')
%% Things to be kept in mind
% datafile name : starlink_satellite_data_id_date,
% spire_satellite_data_id_date
%% user inputs
year_data = '2022';
month_data = '02';
day_data = {'02', '03', '04', '05'};
flag_rho = 'JB08';  
data_pod = [];    % _ucar or []
hasdm_models = [{'2022_HASDM_400-475KM'}, {'2022_HASDM_500-575KM'}, {'2022_HASDM_600-675KM'}];
erp_ceres_datafile = 'CERES_EBAF_Ed4.2_Subset_202202-202202.nc';
case_run = 'EDR';
del_T = 1;
%% Hasdm initialize
td = datetime(str2double(year_data),str2double(month_data),str2double(day_data));
doy =  day(td, 'dayofyear');  % 302; %                     %% day of year (nrlmsise-00)

% center a few days around doy
n_days = range(doy)+2;
doy_hasdm = min(doy)-1;
hasdm_mat = hasdm_initialize(n_days, doy_hasdm, hasdm_models);   % hasdm sub-matrix
%% Point to the data directory and download the data
date_curr = strcat(year_data,'_',month_data,'_',day_data{1});
% dir_data = (strcat('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/main_working_folder_data/spire/2022_02_05', data_pod));
dir_data = (strcat('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods/Simultaneous_estimation', data_pod));
dir_struct = dir(dir_data);
dir_name = {dir_struct.name};
dir_name = dir_name(contains(dir_name, 'satellite_data_'));
sat_ID_mat = extractBetween(dir_name, strcat('data', data_pod, '_'),  strcat('_', year_data)); % {'FM103'}; %
dataset_sat_mat = extractBefore(dir_name, '_satellite');

for ii = 1:numel(sat_ID_mat)
    ii
    sat_ID = sat_ID_mat{ii};
    dataset_sat = dataset_sat_mat{ii};
    run_ODcode(dir_data, sat_ID, date_curr, flag_rho, dataset_sat, hasdm_mat, erp_ceres_datafile, data_pod, case_run, del_T);
end
toc

function [] = run_ODcode(dir_data,sat_ID, date_curr, flag_rho, dataset_sat, hasdm_mat, erp_ceres_datafile, data_pod, case_run, del_T)
Main_simulStudy
end