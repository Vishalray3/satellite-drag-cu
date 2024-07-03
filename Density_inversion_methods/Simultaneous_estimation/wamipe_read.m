%% wam-ipe read and store in mat file
clc
clearvars
path_name = '/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods/data/wamipe_data';
addpath(path_name)
filename = 'wfs.t00z.gsm05.20230415_210500';
rho_wam_exp = ncread(strcat(filename,'.nc'), 'den');
vec_alt = ncread(strcat(filename,'.nc'), 'hlevs');
vec_lat = ncread(strcat(filename,'.nc'), 'lat');
vec_long = ncread(strcat(filename,'.nc'), 'lon');
%% remove non-zero values
n_lat = numel(vec_lat);
vec_lat = vec_lat(2:n_lat-1);
rho_wam = log(rho_wam_exp(:,2:n_lat-1,:,:));
% rho_wam = rho_wam(ind_lon, ind_lat, ind_alt, ind_tim);
% vec_alt = vec_alt(ind_alt);
% vec_lat = vec_lat(ind_lat);
% vec_lon = vec_lon(ind_lon);

%%
% save(strcat(path_name,'/',filename),'vec_alt','vec_lat','vec_long','rho_wam')
