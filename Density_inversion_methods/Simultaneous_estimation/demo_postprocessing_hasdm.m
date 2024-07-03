%% Post processing
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/main_working_folder_data')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods/data/HASDM_data')
clc
clearvars
alt_ml = [];
for ii = 1
    data_struct_jb08 = importdata(strcat('KayhanDensities_JB082022_02_0',num2str(ii),'.csv'));
    data_struct_wamipe = importdata(strcat('KayhanDensities_WAMIPE2022_02_0',num2str(ii),'.csv'));
    data_struct_msis = importdata(strcat('KayhanDensities_MSIS002022_02_0',num2str(ii),'.csv'));
    data_struct2 = importdata('Processed_Densities_02-05-2022.csv');
    
    Altitudejb08 = [data_struct_jb08.data(:,6)];
    Latitudejb08 = [data_struct_jb08.data(:,4)];
    Longitudejb08 = [data_struct_jb08.data(:,5)];
    Estimated_jb08 = [data_struct_jb08.data(:,2)];
    JB1 = [data_struct_jb08.data(:,7)];
    time_jd_jb08 = [data_struct_jb08.data(:,3)];
    sat_id1_jb08 = [data_struct_jb08.data(:,1)];
    
    Altitudewamipe = [data_struct_wamipe.data(:,6)];
    Latitudewamipe = [data_struct_wamipe.data(:,4)];
    Longitudewamipe = [data_struct_wamipe.data(:,5)];
    Estimated_wamipe = [data_struct_wamipe.data(:,2)];
    wamipe = [data_struct_wamipe.data(:,7)];
    time_jd_wamipe = [data_struct_wamipe.data(:,3)];
    sat_id1_wamipe = [data_struct_wamipe.data(:,1)];
    
    Altitudemsis = [data_struct_msis.data(:,6)];
    Latitudemsis = [data_struct_msis.data(:,4)];
    Longitudemsis = [data_struct_msis.data(:,5)];
    Estimated_msis = [data_struct_msis.data(:,2)];
    msis = [data_struct_msis.data(:,7)];
    time_jd_msis = [data_struct_msis.data(:,3)];
    sat_id1_msis = [data_struct_msis.data(:,1)];
    
    time_stamp = [data_struct2.textdata(2:end,1)];
    HASDM_ml = [data_struct2.data(:,11)];
    alt_ml = [data_struct2.data(:,4)];
    sat_id2 = [data_struct2.data(:,1)];
end

%% HASDM densities
hasdm_models = [{'2022_HASDM_175-275KM'}, {'2022_HASDM_300-375KM'}, {'2022_HASDM_400-475KM'}, {'2022_HASDM_500-575KM'}, {'2022_HASDM_600-675KM'}];

doy = 32;
n_days = 3;
doy = doy-1;
hasdm_mat = hasdm_initialize(n_days, doy, hasdm_models);   % hasdm sub-matrix


F_hasdm = hasdm_interpolant(n_days, doy, hasdm_mat);
jd_ref = 86400*GREGORIANtoJD_vector(1950,1,0) + 3600*0 + 60*0 + 0;


njd = (time_jd_jb08 - jd_ref)/86400;
long_pos = Longitudejb08;
long_pos(long_pos<0) = long_pos(long_pos<0) + 360;
alt = Altitudejb08;
lat = Latitudejb08;
HASDM_jb08 = exp(F_hasdm(alt, njd,long_pos, lat));

njd = (time_jd_wamipe - jd_ref)/86400;
long_pos = Longitudewamipe;
long_pos(long_pos<0) = long_pos(long_pos<0) + 360;
alt = Altitudewamipe;
lat = Latitudewamipe;
HASDM_wamipe = exp(F_hasdm(alt, njd,long_pos, lat));

njd = (time_jd_msis - jd_ref)/86400;
long_pos = Longitudemsis;
long_pos(long_pos<0) = long_pos(long_pos<0) + 360;
alt = Altitudemsis;
lat = Latitudemsis;
HASDM_msis = exp(F_hasdm(alt, njd,long_pos, lat));


%%
% jd_hml1 = 86400*juliandate(datetime(time_stamp(1:115536),'InputFormat','yyyy-MM-dd HH:mm:ss+00:00'));
% jd_hml2 = 86400*juliandate(datetime(time_stamp(115537:end),'InputFormat','yyyy-MM-dd HH:mm:ss.S+00:00'));
% jd_hml = [jd_hml1;jd_hml2];
% 
% sat_id_uni = unique(sat_id1_jb08);
% for kk = 1:numel(sat_id_uni)
%     ind_id_hml = find(sat_id2 == sat_id_uni(kk));
%     HASDM_trunc = HASDM_ml(ind_id_hml);
%     alt_trunc = alt_ml(ind_id_hml);
%     jd_trunc = jd_hml(ind_id_hml);
%     
%     ind_id = find(sat_id1_jb08 == sat_id_uni(kk));
%     for ii = 1:numel(ind_id)
%         [~,ind_trunc] = min(abs(time_jd_jb08(ind_id(ii))-jd_hml(ind_id_hml)));
%         HASDM_jb08(ind_id(ii),1) = HASDM_trunc(ind_trunc);
%         alt_new_jb08(ind_id(ii),1) = alt_trunc(ind_trunc);
%         jd_new_jb08(ind_id(ii),1) = jd_trunc(ind_trunc);
%     end
%    
%     ind_id = find(sat_id1_wamipe == sat_id_uni(kk));
%     for ii = 1:numel(ind_id)
%         [~,ind_trunc] = min(abs(time_jd_wamipe(ind_id(ii))-jd_hml(ind_id_hml)));
%         HASDM_wamipe(ind_id(ii),1) = HASDM_trunc(ind_trunc);
%         alt_new_wamipe(ind_id(ii),1) = alt_trunc(ind_trunc);
%         jd_new_wamipe(ind_id(ii),1) = jd_trunc(ind_trunc);
%     end
% 
%     ind_id = find(sat_id1_msis == sat_id_uni(kk));
%     for ii = 1:numel(ind_id)
%         [~,ind_trunc] = min(abs(time_jd_msis(ind_id(ii))-jd_hml(ind_id_hml)));
%         HASDM_msis(ind_id(ii),1) = HASDM_trunc(ind_trunc);
%         alt_new_msis(ind_id(ii),1) = alt_trunc(ind_trunc);
%         jd_new_msis(ind_id(ii),1) = jd_trunc(ind_trunc);
%     end    
% end


date_curr = datetime(time_jd_jb08/86400, 'ConvertFrom', 'juliandate');
edges = [200:25:650];
Yalt_jb08 = discretize(Altitudejb08, edges);
Yalt_wamipe = discretize(Altitudewamipe, edges);
Yalt_msis = discretize(Altitudemsis, edges);


indices = unique(Yalt_jb08);

for ii = 1:numel(indices)
    vec_ind_jb08 = find(Yalt_jb08 == indices(ii));
    rho_hasdm = HASDM_jb08(vec_ind_jb08);
    rho_est_jb08 = Estimated_jb08(vec_ind_jb08);
    rho_jb08 = JB1(vec_ind_jb08);
    error_est_jb08(ii) = (mean(rho_est_jb08)-mean(rho_hasdm))/mean(rho_hasdm)*100;
    error_jb08(ii) = (mean(rho_jb08)-mean(rho_hasdm))/mean(rho_hasdm)*100;
    errseries_est = (rho_est_jb08 - rho_hasdm)./rho_hasdm*100;
    errseries_jb08 = (rho_jb08 - rho_hasdm)./rho_hasdm*100;
    mean_est_jb08(ii) = mean(abs(errseries_est));
    rms_est_jb08(ii)  = rms(errseries_est);
    mean_jb08(ii) = mean(abs(errseries_jb08));
    rms_jb08(ii)  = rms(errseries_jb08);
    
    vec_ind_wamipe = find(Yalt_wamipe == indices(ii));
    rho_hasdm = HASDM_wamipe(vec_ind_wamipe);
    rho_est_wamipe = Estimated_wamipe(vec_ind_wamipe);
    rho_wamipe = wamipe(vec_ind_wamipe);
    error_est_wamipe(ii) = (mean(rho_est_wamipe)-mean(rho_hasdm))/mean(rho_hasdm)*100;
    error_wamipe(ii) = (mean(rho_wamipe)-mean(rho_hasdm))/mean(rho_hasdm)*100;
    errseries_est = (rho_est_wamipe - rho_hasdm)./rho_hasdm*100;
    errseries_wamipe = (rho_wamipe - rho_hasdm)./rho_hasdm*100;
    mean_est_wamipe(ii) = mean(abs(errseries_est));
    rms_est_wamipe(ii)  = rms(errseries_est);
    mean_wamipe(ii) = mean(abs(errseries_wamipe));
    rms_wamipe(ii)  = rms(errseries_wamipe);
    
    vec_ind_msis = find(Yalt_msis == indices(ii));
    rho_hasdm = HASDM_msis(vec_ind_msis);
    rho_est_msis = Estimated_msis(vec_ind_msis);
    rho_msis = msis(vec_ind_msis);
    error_est_msis(ii) = (mean(rho_est_msis)-mean(rho_hasdm))/mean(rho_hasdm)*100;
    error_msis(ii) = (mean(rho_msis)-mean(rho_hasdm))/mean(rho_hasdm)*100;
    errseries_est = (rho_est_msis - rho_hasdm)./rho_hasdm*100;
    errseries_msis = (rho_msis - rho_hasdm)./rho_hasdm*100;
    mean_est_msis(ii) = mean(abs(errseries_est));
    rms_est_msis(ii)  = rms(errseries_est);
    mean_msis(ii) = mean(abs(errseries_msis));
    rms_msis(ii)  = rms(errseries_msis);
    
    
    
end
error_mat = [error_est_jb08; error_jb08;error_est_msis; error_msis;error_est_wamipe; error_wamipe];
mean_mat = [mean_est_jb08; mean_jb08;mean_est_msis; mean_msis;mean_est_wamipe; mean_wamipe];
rms_mat = [rms_est_jb08; rms_jb08;rms_est_msis; rms_msis;rms_est_wamipe; rms_wamipe];


%%
Yalt_mean = [225:25:350, 475:25:650];
figure(1)
b = bar(Yalt_mean, abs(error_mat));
xlabel('Altitude bins (km)')
ylabel('Error (%)')
title('Average Density Difference w.r.t HASDM')
legend('Corrected JB2008','JB2008','Corrected MSIS00','MSIS00','Corrected WAMIPE','WAMIPE')
set(gca, 'FontSize', 16)
grid on
set(gca,'YScale','log')
b(1).FaceColor = [1 0 0];
b(1).LineStyle = "-";
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(2).LineStyle = ":";
b(3).FaceColor = [0 1 1];
b(4).FaceColor = [0.3010 0.7450 0.9330];
b(4).LineStyle = ":";
b(5).FaceColor = [0 1 0];
b(6).FaceColor = [0.4660 0.6740 0.1880];
b(6).LineStyle = ":";
%%
figure(2)
subplot(2,1,1)
b = bar(Yalt_mean, mean_mat);
ylabel('Mean Error (%)')
title('Density difference w.r.t HASDM')
yticks([1, 10, 100])
% legend('Corrected JB2008','JB2008','Corrected MSIS00','MSIS00','Corrected WAMIPE','WAMIPE')
set(gca, 'FontSize', 16)
grid on
set(gca,'YScale','log')
b(1).FaceColor = [1 0 0];
b(1).LineStyle = "-";
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(2).LineStyle = ":";
b(3).FaceColor = [0 1 1];
b(4).FaceColor = [0.3010 0.7450 0.9330];
b(4).LineStyle = ":";
b(5).FaceColor = [0 1 0];
b(6).FaceColor = [0.4660 0.6740 0.1880];
b(6).LineStyle = ":";

subplot(2,1,2)
b = bar(Yalt_mean, rms_mat);
xlabel('Altitude bins (km)')
ylabel('RMS Error (%)')
legend('Corrected JB2008','JB2008','Corrected MSIS00','MSIS00','Corrected WAMIPE','WAMIPE')
% legend('Estimate','JB2008')
set(gca, 'FontSize', 16)
yticks([1, 10, 100])
grid on
set(gca,'YScale','log')
b(1).FaceColor = [1 0 0];
b(1).LineStyle = "-";
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(2).LineStyle = ":";
b(3).FaceColor = [0 1 1];
b(4).FaceColor = [0.3010 0.7450 0.9330];
b(4).LineStyle = ":";
b(5).FaceColor = [0 1 0];
b(6).FaceColor = [0.4660 0.6740 0.1880];
b(6).LineStyle = ":";
%%
% figure(3)
% plot(Estimated_jb08)
% hold on
% plot( HASDM)
% plot(JB1)
% % xlabel('Altitude bins (km)')
% ylabel('Density (kg/m3)')
% title('Estimated density')
% legend('Estimate','HASDM','JB2008')
% set(gca, 'FontSize', 16)
% grid on