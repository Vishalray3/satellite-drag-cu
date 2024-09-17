%% Live demo script
clc
clearvars
% tic
restoredefaultpath
linux_os = 1;
data_pod = [];    % _ucar or []
parallel_flag = 1;
[parent_directory, dir_data, output_dir] = data_paths(linux_os);
%% Things to be kept in mind
% datafile name : starlink_satellite_data_id_date,
% spire_satellite_data_id_date

%% Furnishing spice kernels to use the spice functions for each parallel thread
if parallel_flag == 1
    c = parcluster('local'); % build the 'local' cluster object
    numWorkers = c.NumWorkers;
    parfor i = 1:numWorkers
        % Add paths and furnish kernels
        cspice_furnsh('de430.bsp')                       %% Planetary ephemerides
        cspice_furnsh('naif0012.tls')                    %% time leap seconds
        cspice_furnsh('earth_assoc_itrf93.tf')           %% Earth ITRF frame
        cspice_furnsh('pck00010.tpc')
        cspice_furnsh('gm_de431.tpc')                    %% GM values
        cspice_furnsh('earth_000101_210629_210407.bpc')
    end
else
    cspice_furnsh('de430.bsp')                       %% Planetary ephemerides
    cspice_furnsh('naif0012.tls')                    %% time leap seconds
    cspice_furnsh('earth_assoc_itrf93.tf')           %% Earth ITRF frame
    cspice_furnsh('pck00010.tpc')
    cspice_furnsh('gm_de431.tpc')                    %% GM values
    cspice_furnsh('earth_000101_210629_210407.bpc')
end
%% user inputs
year_data = 2022;
month_sat_array = [1:12];
hasdm_models = [{'2022_HASDM_400-475KM'}, {'2022_HASDM_500-575KM'}, {'2022_HASDM_600-675KM'}];
flag_rho = 'MSIS20';  
erp_ceres_datafile = 'CERES_EBAF_Ed4.2_Subset_202202-202202.nc';
case_run = 'EDR';
del_T = 1;
sat_ids_skip = {'FM088', 'FM084', 'FM085', 'FM087', 'FM088', 'FM099', 'FM100', 'FM102', 'FM103', 'FM104'};
%% Hasdm initialize
den_mat_list = cell(1, numel(hasdm_models));

for ii = 1:numel(hasdm_models)
    load(hasdm_models{ii}, 'den_mat')
    den_mat_list{ii} = den_mat;
end

%% Get the satelllite IDs
dir_struct = dir(dir_data);
dir_name = {dir_struct.name};
data_pattern = "spire_satellite_data_FM" + digitsPattern(3) + '_' + digitsPattern(4) + '_' + digitsPattern(2) + '_' + digitsPattern(2) + '.mat';

dir_name = dir_name(matches(dir_name, data_pattern));
sat_ID_mat = extractBetween(dir_name, strcat('data', data_pod, '_'),  strcat('_', sprintf('%d', year_data))); %  %
dataset_sat_mat = extractBefore(dir_name, '_satellite');
sat_ID_mat(ismember(sat_ID_mat, sat_ids_skip)) = [];
% sat_ID_mat = {'FM102', 'FM103'};
%% Main loop
dataset_sat = dataset_sat_mat{1};
parfor ii = 1:numel(sat_ID_mat)
    sat_ID = sat_ID_mat{ii};
    sat_ID
    run_parfor_loop(dir_data, month_sat_array, year_data, den_mat_list, sat_ID, flag_rho, dataset_sat, erp_ceres_datafile, data_pod, case_run, del_T, linux_os, output_dir)
end 
% toc


function run_parfor_loop(dir_data, month_sat_array, year_data, den_mat_list, sat_ID, flag_rho, dataset_sat, erp_ceres_datafile, data_pod, case_run, del_T, linux_os, output_dir)
    for month_data = month_sat_array
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
            try
                run_ODcode(dir_data, sat_ID, date_curr, flag_rho, dataset_sat, hasdm_mat, erp_ceres_datafile, data_pod, case_run, del_T, linux_os, output_dir);
            catch error_loop
%                 msg_txt = getReport(error_loop);
%                 fprintf(1,'There was an error for satellite %s on day %u and month %u! The message was:%s',sat_ID, day_data, month_data, msg_txt); 
                fprintf(1,'There was an error for satellite %s on day %u and month %u! The message was:\n%s',sat_ID, day_data, month_data, error_loop.message); 
            end
        end
    end 
end

function [] = run_ODcode(dir_data,sat_ID, date_curr, flag_rho, dataset_sat, hasdm_mat, erp_ceres_datafile, data_pod, case_run, del_T, linux_os, output_dir)
Main_simulStudy
end