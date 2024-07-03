close all
clc
clearvars
%%
year_data = '2022';
month_data = '02';
day_data = '05';
flag_rho = 'JB08';  
data_pod = '_ucar';    % '_ucar' or []
dir_data = (strcat('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods/data/spire/2022_02_05', data_pod));
dir_struct = dir(dir_data);
dir_name = {dir_struct.name};
dir_name = dir_name(contains(dir_name, 'satellite_data_'));
% sat_ID_mat = extractBetween(dir_name, strcat('data', data_pod, '_'),  strcat('_', year_data)); % {'FM103'}; %
sat_ID_mat = {'FM099','FM100','FM102','FM103','FM104','FM106','FM117','FM118','FM119','FM120','FM122','FM124','FM125','FM128','FM129'};
date_curr = strcat(year_data,'_',month_data,'_',day_data);

mu_e = 398600435436096;
%%
rho_true_all = [];
rho_est_all = [];
rho_nom_all = [];
alt_all = [];

for ii = 1: numel(sat_ID_mat)
    sat_ID = sat_ID_mat{ii};
    load(strcat(dir_data,'/spireResults_', flag_rho, '_',sat_ID,'_', date_curr, '_new'))
    %% Orbit period
    r_circ = a_sma;
    v_circ = sqrt(mu_e/r_circ);
    n_mean = v_circ/r_circ;
    T_orb = 2*pi/n_mean;
    %%
    arclen = T_orb;
    
    rho_true = rho_hasdm_trunc;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Need to investigate the Cd variation due to area/attitude - makes all
    % the difference
    rho_est = Cd_est(ind_vec).*rho_est_trunc/0.31;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rho_nom = rho_nom_trunc;
    
    time_diff = diff(time_prop_trunc);
    time_ind = find(time_diff > 100);
    time_ind_all = [1, time_ind, numel(time_prop_trunc)];
    
    
    time_step = time_prop_trunc(2) - time_prop_trunc(1);
    
    [rho_true_results, rho_nom_results, rho_est_results, alt_results] = den_avg_spire(rho_true, rho_est, rho_nom, arclen, time_step, time_ind_all, alt_trunc);
    % remove outliers
    [~, ind_out] = rmoutliers(rho_est_results./rho_nom_results);
    
    ind_retain = ~ind_out;
    rho_true_all = [rho_true_all, rho_true_results(ind_retain)];
    rho_est_all = [rho_est_all, rho_est_results(ind_retain)];
    rho_nom_all = [rho_nom_all, rho_nom_results(ind_retain)];
    alt_all = [alt_all, alt_results(ind_retain)];
    
end
alt_all = alt_all*1e-3;
edges = [400:50:650];
Yalt_spire = discretize(alt_all, edges);
indices = unique(Yalt_spire);

for ii = 1:numel(indices)
    vec_ind_alt = find(Yalt_spire == indices(ii));
    rho_hasdm = rho_true_all(vec_ind_alt);
    rho_est = rho_est_all(vec_ind_alt);
    rho_nom = rho_nom_all(vec_ind_alt);
    error_est_jb08(ii) = (mean(rho_est)-mean(rho_hasdm))/mean(rho_hasdm)*100;
    error_jb08(ii) = (mean(rho_nom)-mean(rho_hasdm))/mean(rho_hasdm)*100;
    errseries_est = (rho_est - rho_hasdm)./rho_hasdm*100;
    errseries_jb08 = (rho_nom - rho_hasdm)./rho_hasdm*100;
    mean_est_jb08(ii) = mean(abs(errseries_est));
    rms_est_jb08(ii)  = rms(errseries_est);
    mean_jb08(ii) = mean(abs(errseries_jb08));
    rms_jb08(ii)  = rms(errseries_jb08);
    num_data(ii) = numel(vec_ind_alt);
      
end
error_mat = [error_est_jb08; error_jb08];
mean_mat = [mean_est_jb08; mean_jb08];
rms_mat = [rms_est_jb08; rms_jb08];
%% PLots


Yalt_mean = edges(indices);
figure(1)
b = bar(Yalt_mean, abs(error_mat));
xlabel('Altitude bins (km)')
ylabel('Error (%)')
title('Average Density Difference w.r.t HASDM (UCAR)')
legend('Corrected JB2008','JB2008')
set(gca, 'FontSize', 16)
grid on



% 
figure(2)
subplot(2,1,1)
b = bar(Yalt_mean, mean_mat);
ylabel('Mean Error (%)')
title('Density difference w.r.t HASDM (UCAR)')
% legend('Corrected JB2008','JB2008','Corrected MSIS00','MSIS00','Corrected WAMIPE','WAMIPE')
set(gca, 'FontSize', 16)
grid on

subplot(2,1,2)
b = bar(Yalt_mean, rms_mat);
xlabel('Altitude bins (km)')
ylabel('RMS Error (%)')
legend('Corrected JB2008','JB2008')
% legend('Estimate','JB2008')
set(gca, 'FontSize', 16)
grid on

