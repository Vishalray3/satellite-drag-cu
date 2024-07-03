%% Post-processing
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods/data/wamipe_data')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods/data/hasdm_data')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/JB08')
clc
clearvars
%% Inputs
year_data = '2022';
month_data = '02';
day_data_all = {'05'};
n_days = 1;
flag_rho = 'WAMIPE';
wam_filename = 'WAM_den_20220205';
data_mat = [];

%% Extract datafiles
for ii = 1:numel(day_data_all)
    day_data = day_data_all{ii};
date_curr = strcat(year_data,'_',month_data,'_',day_data);
dir_data = strcat('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods/data/demo');%,date_curr);
addpath(dir_data)
dir_struct = dir(dir_data);
dir_name = {dir_struct.name};
% dir_name = dir_name(contains(dir_name, 'spireResults_JB08'));
dir_name = dir_name(contains(dir_name, strcat('Results_',flag_rho)));
% sat_ID_mat = extractBetween(dir_name, 'FM', '_2019');
sat_ID_mat = extractBetween(dir_name, strcat(flag_rho,'_'), strcat('_',year_data));
sat_ID_mat = erase(sat_ID_mat, 'FM');
sat_id_num = str2double(sat_ID_mat);
%% Read WAM-IPE data (if needed)
% F_wamipe = wamipe_interpolant(wam_filename, str2double(year_data), str2double(month_data), str2double(day_data));
% jd_ref = 86400*GREGORIANtoJD_vector(1950,1,0) + 3600*0 + 60*0 + 0;
%% JB08 Densities
load SOLFSMY
load DTCFILE
earth_model = wgs84Ellipsoid;                  %% shape model for geodetic to geocentric latitude
flattening = earth_model.Flattening;           %% flattening factor
ind_sol = find(SOLdata(1,:) == str2double(year_data),1);        %% index to point towards data for doy
ind_mag = find(DTCdata(1,:) == str2double(year_data),1);                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HASDM densities
td = datetime(str2double(year_data), str2double(month_data), str2double(day_data));
doy =  day(td, 'dayofyear');  % 302; %  
F_hasdm = hasdm_interpolant(n_days, doy);
%% Write data
for ii = 1:numel(dir_name)
    ii
    sat_ID = sat_id_num(ii);
    load(dir_name{ii})
    if ~isempty(rho_est_trunc)
        sat_ID = sat_ID*ones(1,numel(rho_est_trunc));
%         rho_hasdm_trunc = zeros(1, numel(rho_est_trunc));
%         njd = (jdutc_trunc - jd_ref)/86400;
        long_trunc_pos = long_trunc;
        long_trunc_pos(long_trunc_pos<0) = long_trunc_pos(long_trunc_pos<0) + 360;
        rho_wamipe_trunc = zeros(1, numel(rho_est_trunc)); %exp(F_wamipe(long_trunc_pos, lat_trunc,alt_trunc/1e3, njd));
        rho_hasdm_trunc = zeros(1, numel(rho_est_trunc)); %exp(F_hasdm(alt_trunc/1e3, njd,long_trunc_pos, lat_trunc));
%         rho_jb08_trunc = rho_nom_trunc; %density_jb08(epoch, doy,yyyy(1), time_prop_trunc, alt_trunc, lat_trunc,flattening, X_state_trunc, sun_pos_trunc,...
%             SOLdata, DTCdata, ind_sol, ind_mag);
        data_id = [sat_ID', rho_est_trunc', jdutc_trunc',lat_trunc', long_trunc', alt_trunc'/1e3, rho_nom_trunc', rho_hasdm_trunc', rho_wamipe_trunc'];
        data_mat = [data_mat;data_id];
    end
end
end
data_tab = array2table(data_mat);
data_tab.Properties.VariableNames = {'Sat ID','Estimated Density (kg/m3)','UTC Julian date (sec)','Latitude (deg)','Longitude (deg)',...
    'Altitude (km)','JB2008', 'HASDM','WAM-IPE'};
date_curr = strcat(year_data,'_',month_data,'_',day_data);
writetable(data_tab, strcat('KayhanDensities_',flag_rho,date_curr,'.csv'))

