%% Parallel density estimation from spire data
% Runs density estimation code for multiple satellites independenty

% Have a User input file specified - dates, satellite ID, run
% configuration. 

% Read Spire data and form spire datasets
%% Comments
% what happens during data gaps? should we not use rho_est and use rho_nom
% instead?
% Include ERP
% include center of mass shift

%%
clc
clearvars
tic 
%% user inputs
date_curr = '2022_02_03';
flag_rho = 'MSIS00';   
dataset_sat = 'starlink';
%% Point to the data directory and download the data
dir_data = strcat('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods/data/starlink/',date_curr);
dir_struct = dir(dir_data);
dir_name = {dir_struct.name};
dir_name = dir_name(contains(dir_name, 'starlink_data_'));
sat_ID_mat = extractBetween(dir_name, 'data_', '_2022'); % {'FM103'}; %

for ii = 1:numel(sat_ID_mat)
    sat_ID = sat_ID_mat{ii};
    run_ODcode(dir_data, sat_ID, date_curr, flag_rho, dataset_sat);
end
toc

function [] = run_ODcode(dir_data,sat_ID, date, flag_rho, dataset_sat)
Main_simulStudy
end