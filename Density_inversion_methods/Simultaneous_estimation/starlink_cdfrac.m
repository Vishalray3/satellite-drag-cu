%% STARLINK CD calculation
% Cd_new*rho_hasdm = rho_est*Cd_old
clc
clearvars
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods/data/starlink/2022_02_03')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods/data/starlink')
data_table = readtable('starlink3167_ephemeris_full_hasdm.csv');
rho_hasdm_all = table2array(data_table(:,23));
rho_hasdm = rho_hasdm_all(100:699)';

load('starlinkResults_JB08_3167_2022_02_03')

Cd_new_jb08 = rho_est.*Cd_est./rho_hasdm;

f_new_jb08 = (mean(Cd_new_jb08) - mean(Cd_s))./(mean(Cd_ads) - mean(Cd_s));

rho_new = rho_est.*Cd_est./(f_new_jb08*Cd_ads + (1-f_new_jb08)*Cd_s);
err_nom_jb08 = (rho_nom-rho_hasdm)./rho_hasdm*100;
err_new_jb08 = (rho_new-rho_hasdm)./rho_hasdm*100;
 
err_est_jb08 = (rho_est-rho_hasdm)./rho_hasdm*100;
load('starlinkResults_MSIS00_3167_2022_02_03')

Cd_new_msis = rho_est.*Cd_est./rho_hasdm;

f_new_msis = (mean(Cd_new_msis) - mean(Cd_s))./(mean(Cd_ads) - mean(Cd_s));

% coming out to be 0.28-0.29 for both MSIS and jb2008 and two different
% arcs.

rho_new = rho_est.*Cd_est./(f_new_msis*Cd_ads + (1-f_new_msis)*Cd_s);
err_nom_msis = (rho_nom-rho_hasdm)./rho_hasdm*100;
err_new_msis = (rho_new-rho_hasdm)./rho_hasdm*100;
 
err_est_msis = (rho_est-rho_hasdm)./rho_hasdm*100;

figure(1)
plot(abs(err_nom_msis))
hold on
plot(abs(err_est_msis))
plot(abs(err_new_msis))

figure(2)
plot(abs(err_nom_jb08))
hold on
plot(abs(err_est_jb08))
plot(abs(err_new_jb08))