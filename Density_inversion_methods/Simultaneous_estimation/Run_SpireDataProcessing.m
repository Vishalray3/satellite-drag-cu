%% Live demo script
clc
clearvars
% tic
restoredefaultpath
linux_os = 1;
data_pod = [];    % _ucar or []

if linux_os == 1
    parent_directory = '/home/vira0155';
    dir_data = '/media/faraday/DATA/thermospheric/spire_data/2022/spire_matlab';
    output_dir = fullfile(dir_data, 'results');
else
    parent_directory = '/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder';
    dir_data = fullfile(parent_directory, 'satellite-drag-cu/Density_inversion_methods/Simultaneous_estimation');
    output_dir = dir_data;
end

addpath(fullfile(parent_directory, 'satellite-drag-cu/Density_inversion_methods/data/HASDM_data'))
%% Things to be kept in mind
% datafile name : starlink_satellite_data_id_date,
% spire_satellite_data_id_date
%% user inputs
year_data = 2022;
month_sat_array = [1:12];
hasdm_models = [{'2022_HASDM_400-475KM'}, {'2022_HASDM_500-575KM'}, {'2022_HASDM_600-675KM'}];
flag_rho = 'MSIS00';  
erp_ceres_datafile = 'CERES_EBAF_Ed4.2_Subset_202202-202202.nc';
case_run = 'EDR';
del_T = 1;
%% Hasdm initialize
den_mat_list = cell(1, numel(hasdm_models));

for ii = 1:numel(hasdm_models)
    load(hasdm_models{ii}, 'den_mat')
    den_mat_list{ii} = den_mat;
end

%% Main loop
parfor month_data = month_sat_array
    month_data
    day_sat_array = [1:eomday(year_data, month_data)];
    
    for day_data = day_sat_array
        day_data
        td = datetime(year_data, month_data, day_data);
        doy =  day(td, 'dayofyear');  % 302; %                     %% day of year (nrlmsise-00)

        % center a few days around doy
        n_days = 3;
        doy_hasdm = doy-1;
        hasdm_mat = hasdm_initialize(n_days, doy_hasdm, den_mat_list);   % hasdm sub-matrix
        %% Point to the data directory and download the data
        date_curr = strcat(sprintf('%d', year_data),'_', sprintf('%02d', month_data),'_', sprintf('%02d',day_data));
        dir_struct = dir(dir_data);
        dir_name = {dir_struct.name};
        data_pattern = "spire_satellite_data_FM" + digitsPattern(3) + '_' + digitsPattern(4) + '_' + digitsPattern(2) + '_' + digitsPattern(2) + '.mat';

        dir_name = dir_name(matches(dir_name, data_pattern));
        sat_ID_mat = {'FM100'}; %extractBetween(dir_name, strcat('data', data_pod, '_'),  strcat('_', sprintf('%d', year_data))); %  %
        dataset_sat_mat = extractBefore(dir_name, '_satellite');

        for ii = 1:numel(sat_ID_mat)
            sat_ID = sat_ID_mat{ii};
            dataset_sat = dataset_sat_mat{ii};
            try
                run_ODcode(dir_data, sat_ID, date_curr, flag_rho, dataset_sat, hasdm_mat, erp_ceres_datafile, data_pod, case_run, del_T, linux_os, output_dir);
            catch error_loop
                msg_txt = getReport(error_loop);
                fprintf(1,'There was an error for satellite %s on day %u and month %u! The message was:%s',sat_ID, day_data, month_data, msg_txt); 
%                 fprintf(1,'There was an error for satellite %s on day %u and month %u! The message was:\n%s',sat_ID, day_data, month_data, error_loop.message); 
            end
        end
    end
end 
% toc

function [] = run_ODcode(dir_data,sat_ID, date_curr, flag_rho, dataset_sat, hasdm_mat, erp_ceres_datafile, data_pod, case_run, del_T, linux_os, output_dir)
Main_simulStudy
end