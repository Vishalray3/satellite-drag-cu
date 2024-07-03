%% Post-processing
write_to_api(sprintf("'Writing processed data to file'"),flag_disp)
%% Inputs
year_data = '2022';
month_data = '02';
day_data = '05';
wam_filename = 'WAM_den_20220205';
%% Extract datafiles 
date_curr = strcat(year_data,'_',month_data,'_',day_data);
addpath(dir_data)
dir_struct = dir(dir_data);
dir_name = {dir_struct.name};
% dir_name = dir_name(contains(dir_name, 'spireResults_JB08'));
dir_name = dir_name(contains(dir_name, 'ResultsRho_JB08'));
% sat_ID_mat = extractBetween(dir_name, 'FM', '_2019');
sat_ID_mat = extractBetween(dir_name, 'JB08_', '_2022');
sat_ID_mat = erase(sat_ID_mat, 'FM');
sat_id_num = str2double(sat_ID_mat);
data_mat = [];
%% Write data
for ii = 1:numel(dir_name)
    sat_ID = sat_id_num(ii); 
    load(dir_name{ii})
    if ~isempty(jdutc_new)
    sat_ID = sat_ID*ones(1,numel(jdutc_new));
    rho_hasdm_trunc = zeros(1, numel(jdutc_new));
    rho_jb08_trunc = zeros(1, numel(jdutc_new));
    rho_wamipe = zeros(1, numel(jdutc_new));
    rho_est_trunc = zeros(1, numel(jdutc_new));
    data_id = [sat_ID', rho_est_trunc', jdutc_new',latitude', longitude', alt_calc'/1e3, rho_jb08_trunc', rho_hasdm_trunc', rho_wamipe'];
    data_mat = [data_mat;data_id];
    end
end
data_tab = array2table(data_mat);
data_tab.Properties.VariableNames = {'Sat ID','Estimated Density (kg/m3)','UTC Julian date (sec)','Latitude (deg)','Longitude (deg)',...
    'Altitude (km)','JB2008', 'HASDM','WAM-IPE'};
write_to_api(sprintf("'Securely transferring file to SET'"),flag_disp)
% writetable(data_tab, '/var/sftp/spire/KayhanDensities_2022_02_05_full.csv')
writetable(data_tab, '/var/sftp/spire/KayhanDensities_2022_02_05_full.csv')
write_to_api(sprintf("'Kayhan Space hand-off complete. (Upto a delay of two minutes for file acquisition by SET)'"),flag_disp)

