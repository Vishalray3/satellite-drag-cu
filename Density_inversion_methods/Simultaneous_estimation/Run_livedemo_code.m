%% Live demo script
clc
clearvars

addpath('/home/azureuser/mice/mice/src/mice/')
addpath('/home/azureuser/mice/mice/lib/' )
addpath('/home/azureuser/main_hasdm_demo/JB08')
addpath('/home/azureuser/main_hasdm_demo/Density_inversion_methods')
addpath('/home/azureuser/main_hasdm_demo/Density_inversion_methods/data/ancillary_data')
addpath('/home/azureuser/main_hasdm_demo/Density_inversion_methods/data/wamipe_data')
dir_data = '/home/azureuser/main_hasdm_demo/Density_inversion_methods/data/demo';
% addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/DataTools/mice/mice/src/mice/')
% addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/DataTools/mice/mice/lib/' )
% addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/JB08')
% addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods')
% addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods/data/ancillary_data')
% addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods/data/wamipe_data')
% dir_data = '/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods/data/demo';
%% Things to be kept in mind
% parpool('local', 16)
%% user inputs
year_data = '2022';
month_data = '02';
day_data = '05';
flag_rho = 'JB08';  
flag_disp = 1;       % 1: anything should be written to API
%% Point to the data directory and download the data
date_curr = strcat(year_data,'_',month_data,'_',day_data);
dir_struct = dir(dir_data);
dir_name = {dir_struct.name};
dir_name = dir_name(contains(dir_name, 'satellite_data_'));
sat_ID_mat = extractBetween(dir_name, 'data_', strcat('_', year_data)); % {'FM103'}; %
dataset_sat_mat = extractBefore(dir_name, '_satellite');

write_to_api(sprintf("'Downloading Spire and Starlink data files for %s'",[year_data,'-',month_data,'-',day_data]),flag_disp)
pause(8.0)
write_to_api(sprintf("'Reading Spire and Starlink data files into formats required by processing toolkit'"),flag_disp)
pause(4.0)
for ii = 2%:numel(sat_ID_mat)
    ii
    sat_ID = sat_ID_mat{ii};
    dataset_sat = dataset_sat_mat{ii};
    run_ODcode(dir_data, sat_ID, date_curr, flag_rho, dataset_sat, flag_disp);
end
write_data_demo


function [] = run_ODcode(dir_data,sat_ID, date, flag_rho, dataset_sat, flag_disp)
Main_simulStudy
end