%% Starlink postprocessing
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Density inversion methods')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Density inversion methods/Simultaneous estimation/Starlink_group4_7_first_days')
clc
clearvars
%% Read the MSIS2 data
data_table = readtable('starlink3167_ephemeris_full_MSIS2p0.csv');
rho_msis2(:,1) = table2array(data_table(:,8));
rho_msis2_ind = rho_msis2([100:699, 930:2600]);
data_table = readtable('starlink3167_ephemeris_full_hasdm.csv');
rho_hasdm(:,1) = table2array(data_table(:,23));
rho_hasdm_ind = rho_hasdm([100:699, 930:2600]);
data_table = readtable('starlink3167_ephemeris_full.csv');
Cdcalc = table2array(data_table(:,10));
Acalc = table2array(data_table(:,11));
%%
load ('starlink_3167_jb08_preFeb4')
ind_vec1 = index_vec;
t_vec1 = t_vec;
altitude_mat = altitude_ind;
recef_mat = recef_ind;
vecef_mat = vecef_ind;
rho_est1 = rho_est;
rho_nom1 = rho_nom;
jdutc_sec1 = jdutc_sec;
X_eci_mat = X_eci1;
load ('starlink_3167_jb08_postFeb4')
ind_vec1 = [ind_vec1,index_vec];
t_vec1 = [t_vec1;t_vec];
altitude_mat = [altitude_mat;altitude_ind];
recef_mat = [recef_mat;recef_ind];
vecef_mat = [vecef_mat;vecef_ind];
rho_est1 = [rho_est1, rho_est];
rho_nom1 = [rho_nom1, rho_nom];
jdutc_sec1 = [jdutc_sec1, jdutc_sec];
X_eci_mat = [X_eci_mat, X_eci1];

load ('starlink_3167_jb08_preFeb4_lift')
rho_est2 = rho_est;
load ('starlink_3167_jb08_postFeb4_lift')
rho_est2 = [rho_est2, rho_est];

load ('starlink_3167_jb08_preFeb4_quasiCd')
rho_est3 = rho_est;
Cd_est3 = Cd_est;
load ('starlink_3167_jb08_postFeb4_quasiCd')
rho_est3 = [rho_est3, rho_est];
Cd_est3 = [Cd_est3, Cd_est];

for ii = 1:numel(jdutc_sec1)
    X_lla = ecef2lla(recef_mat(ii,:));
    latitude(ii) = X_lla(1);
    longitude(ii) = X_lla(2);
    if longitude(ii) < 0
        longitude(ii) = longitude(ii) + 3600;
    end
    height_calc(ii) = X_lla(3)/1e3;
    coe(ii,:) = rv2coe_E(X_eci_mat(1:3,ii),X_eci_mat(4:6,ii),mu_e);
    Torb_mat(ii) = 2*pi/sqrt(mu_e/coe(ii,1)^3);
end

hour_ind = hour(t_vec1) + minute(t_vec1)/60 + second(t_vec1)/3600;
local_time = hour_ind + longitude*4/60;
prec_rate = -3/2*Re^2/(coe(1,1)*(1-coe(1,2)^2))^2*1.0826e-3*2*pi/Torb_mat(1)*cosd(coe(1,3));

kk = 1;
for ii = 1:numel(jdutc_sec1)
    jdtime_oneorb = jdutc_sec1(ii) + Torb_mat(ii);
    [~, indt] = min(abs(jdutc_sec1 - jdtime_oneorb));
    if jdtime_oneorb >  jdutc_sec1(end)
        break
    end
    if abs(jdutc_sec1(indt) - jdtime_oneorb) < Torb_mat(ii)/2
        [alt_min, ind_sort] =  min(height_calc(ii:indt));
        if alt_min > 211.5
            continue
        end
        ind_sort = ii + ind_sort-1;
        alt_sort(kk) = height_calc(ind_sort);
        lat_sort(kk) = latitude(ind_sort);
        long_sort(kk) = longitude(ind_sort);
        rho_est_sort(kk) = rho_est1(ind_sort);
        rho_est_sort2(kk) = rho_est2(ind_sort);
        rho_est_sort3(kk) = rho_est3(ind_sort);
        Cd_est_sort3(kk) = Cd_est3(ind_sort);
        rho_nom_sort(kk) = rho_nom1(ind_sort);
        rho_msis2_sort(kk) = rho_msis2_ind(ind_sort);
        rho_hasdm_sort(kk) = rho_hasdm_ind(ind_sort);
        t_vec_sort(kk) = t_vec1(ind_sort);
        local_time_sort(kk) = local_time(ind_sort);
        jd_sort(kk) = jdutc_sec1(ind_sort);
        kk = kk+1;
    end
    
end
jd_diff = diff(jd_sort);
ind_diff = find(jd_diff> 1000);
ind_sort = [ind_diff, numel(jd_diff)];
alt_sort = alt_sort(ind_sort);
lat_sort = lat_sort(ind_sort);
long_sort = long_sort(ind_sort);
rho_est_sort = rho_est_sort(ind_sort);
rho_est_sort2 = rho_est_sort2(ind_sort);
rho_est_sort3 = rho_est_sort3(ind_sort);
Cd_est_sort3 = Cd_est_sort3(ind_sort);
rho_nom_sort = rho_nom_sort(ind_sort);
rho_msis2_sort = rho_msis2_sort(ind_sort);
rho_hasdm_sort = rho_hasdm_sort(ind_sort);
t_vec_sort = t_vec_sort(ind_sort);
local_time_sort = local_time_sort(ind_sort);
jd_sort = jd_sort(ind_sort);
%
% clear year
% data_mat = [year(t_vec1),month(t_vec1),day(t_vec1),hour(t_vec1),minute(t_vec1),second(t_vec1),latitude', longitude', height_calc'];
% data_tab = array2table(data_mat);
% data_tab.Properties.VariableNames = {'Year','Month','Day','Hour(UTC)','Minute(UTC)','Second(UTC)','Latitude (deg)','Longitude (deg)','Altitude (km)'};
% writetable(data_tab, 'starlink3167_ephemeris.csv')


%% plots
figure(1)
subplot(2,1,1)
plot(t_vec1,rho_nom1, 'k.')
hold on
plot(t_vec1,rho_msis2_ind, 'g.')
plot(t_vec1,rho_est1, 'r.','LineWidth', 1)
plot(t_vec1,rho_est1, 'c.','LineWidth', 1)
title('Density (kg/m^3)')
legend('JB2008','MSIS2','Starlink','Starlink (w/ lift)')
% title('Density estimate')
set(gca,'FontSize',16)
grid on
subplot(2,1,2)
plot(t_vec1,height_calc, 'k.')
title('Altitude (km)')
set(gca,'FontSize',16)
grid on

%%
figure(2)
% subplot(2,1,1)
plot(t_vec_sort,rho_nom_sort, 'k.', 'MarkerSize', 15)
hold on
plot(t_vec_sort,rho_msis2_sort, 'g.', 'MarkerSize', 15)
plot(t_vec_sort,rho_hasdm_sort, 'r.', 'MarkerSize', 15)
plot(t_vec_sort,rho_est_sort2, 'c*', 'MarkerSize', 5)
plot(t_vec_sort,rho_est_sort3, 'm*', 'MarkerSize', 5)
title('Density near perigee (kg/m^3)')
legend('JB2008','MSIS2','HASDM','Starlink (diffuse Cd)','Starlink (specular Cd)')
% title('Density estimate')
set(gca,'FontSize',14)
ylim([1e-10 4e-10])
% ax = gca;
% ax.XTickLabelMode = 'manual';
grid on
% subplot(2,1,2)
% plot(t_vec_sort,alt_sort, 'k.', 'MarkerSize', 15)
% title('Altitude (km)')
% set(gca,'FontSize',14)
% grid on

%%
figure(3)
plot(t_vec1,Cdcalc(ind_vec1), 'k.')
hold on
plot(t_vec1,Cd_est3./Acalc(ind_vec1)', 'r.')
title('Drag-coefficient')
legend('Diffuse','Specular')
% title('Density estimate')
set(gca,'FontSize',16)
grid on

figure(4)
plot(t_vec1,Acalc(ind_vec1), 'k.')
title('Cross-sectional area')
% legend('JB2008','MSIS2','Starlink','Starlink (w/ lift)')
% title('Density estimate')
set(gca,'FontSize',16)
grid on

