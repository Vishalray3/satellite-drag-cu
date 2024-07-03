%% Post processing
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods/data/Kayhan densities')
clc
clearvars
Altitudekm = [];
EstimatedDensitykgm3 = [];
JB1 = [];
HASDM = [];
time_jd = [];
for ii = 1:7
data_struct = importdata(strcat('KayhanDensities_2019_12_0',num2str(ii),'.csv'));
Altitudekm = [Altitudekm;data_struct.data(:,6)];
EstimatedDensitykgm3 = [EstimatedDensitykgm3;data_struct.data(:,2)];
JB1 = [JB1;data_struct.data(:,7)];
HASDM = [HASDM;data_struct.data(:,8)];
time_jd = [time_jd;data_struct.data(:,3)];
end


date_curr = datetime(time_jd/86400, 'ConvertFrom', 'juliandate');
edges = [460:10:580];
Yalt = discretize(Altitudekm, edges); 

indices = unique(Yalt);

for ii = 1:numel(indices)
    vec_ind = find(Yalt == indices(ii));
    rho_est = EstimatedDensitykgm3(vec_ind);
    rho_jb08 = JB1(vec_ind);
    rho_hasdm = HASDM(vec_ind);
    error_est(ii) = (mean(rho_est)-mean(rho_hasdm))/mean(rho_hasdm)*100;
    error_jb08(ii) = (mean(rho_jb08)-mean(rho_hasdm))/mean(rho_hasdm)*100;
    errseries_est = (rho_est - rho_hasdm)./rho_hasdm*100;
    errseries_jb08 = (rho_jb08 - rho_hasdm)./rho_hasdm*100;
    mean_est(ii) = mean(abs(errseries_est));
    rms_est(ii)  = rms(errseries_est);
    mean_jb08(ii) = mean(abs(errseries_jb08));
    rms_jb08(ii)  = rms(errseries_jb08);
end
error_mat = [error_est; error_jb08];
mean_mat = [mean_est; mean_jb08];
rms_mat = [rms_est; rms_jb08];


%%
Yalt_mean = [475:10:575];
figure(1)
b = bar(Yalt_mean, error_mat);
xlabel('Altitude bins (km)')
ylabel('Error (%)')
title('Error in daily averaged density w.r.t HASDM')
legend('Estimate','JB2008')
set(gca, 'FontSize', 16)
grid on

figure(2)
subplot(2,1,1)
b = bar(Yalt_mean, mean_mat);
xlabel('Altitude bins (km)')
ylabel('Error (%)')
title('Density error mean')
legend('Estimate','JB2008')
set(gca, 'FontSize', 16)
grid on

subplot(2,1,2)
b = bar(Yalt_mean, rms_mat);
xlabel('Altitude bins (km)')
ylabel('Error (%)')
title('Density error rms')
% legend('Estimate','JB2008')
set(gca, 'FontSize', 16)
grid on

%%
figure(3)
plot(EstimatedDensitykgm3)
hold on
plot( HASDM)
plot(JB1)
% xlabel('Altitude bins (km)')
ylabel('Density (kg/m3)')
title('Estimated density')
legend('Estimate','HASDM','JB2008')
set(gca, 'FontSize', 16)
grid on