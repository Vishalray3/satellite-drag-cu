%% load variables and plot the results
clc
clearvars
% tic
restoredefaultpath
linux_os = 1;
data_pod = [];    % _ucar or []

[parent_directory, dir_data, output_dir] = data_paths(linux_os);
fig_path = fullfile(output_dir, 'figures');

%% User inputs
sat_id = "FM102";
year = "2022";
%%
dir_struct = dir(output_dir);
dir_name = {dir_struct.name};
data_pattern = "results_cdpanel_spire_satellite_data_" + sat_id + '_' + year + '_' + digitsPattern(2) + '_' + digitsPattern(2) + '.mat';

file_names = dir_name(matches(dir_name, data_pattern));

rho_data_all = [];
rho_nom_all = [];
rho_hasdm_all = [];
time_rho_all = [];

for ii=1:numel(file_names)
    load(fullfile(output_dir, file_names{ii}))
    rho_hasdm_all = [rho_hasdm_all, rho_hasdm_eff];
    rho_nom_all = [rho_nom_all, rho_nom_eff];
    rho_data_all = [rho_data_all, rho_eff];
    time_rho_all = [time_rho_all, time_rho];
end

%%
nom_ratio = rho_nom_all./rho_hasdm_all;
data_ratio = rho_data_all./rho_hasdm_all;

nom_mean = nanmean(nom_ratio, 2);
nom_rms = rms(nom_ratio', "omitnan")';
data_mean = nanmean(data_ratio, 2);
data_rms = rms(data_ratio', "omitnan")';


figure(1)
subplot(3,1,1)
plot(time_rho_all(1,:)/3600,rho_hasdm_all(1,:),'k.','LineWidth',1)
hold on
plot(time_rho_all(1,:)/3600,rho_nom_all(1,:),'r.','LineWidth',1)
plot(time_rho_all(1,:)/3600,rho_data_all(1,:),'g.','LineWidth',1)
grid on
title('30-minute arc-length')
legend('HASDM','MSIS00','Spire EDR')
set(gca,'FontSize',14)

subplot(3,1,2)
plot(time_rho_all(2,:)/3600,rho_hasdm_all(2,:),'k.','LineWidth',1)
hold on
plot(time_rho_all(2,:)/3600,rho_nom_all(2,:),'r.','LineWidth',1)
plot(time_rho_all(2,:)/3600,rho_data_all(2,:),'g.','LineWidth',1)
grid on
ylabel('Density ($kg/m^3$)','Interpreter','latex')
title('60-minute arc-length')
set(gca,'FontSize',14)

subplot(3,1,3)
plot(time_rho_all(3,:)/3600,rho_hasdm_all(3,:),'k.','LineWidth',1)
hold on
plot(time_rho_all(3,:)/3600,rho_nom_all(3,:),'r.','LineWidth',1)
plot(time_rho_all(3,:)/3600,rho_data_all(3,:),'g.','LineWidth',1)
grid on
ylabel('Density ($kg/m^3$)','Interpreter','latex')
title('90-minute arc-length')
set(gca,'FontSize',14)
saveas(gcf, fullfile(fig_path, strcat('spire_', year, '_', sat_id, '_density_edr.png')))    
    
figure(2)
subplot(3,1,1)
plot(time_rho_all(1,:)/3600, nom_ratio(1,:),'k.','LineWidth',1)
hold on
plot(time_rho_all(1,:)/3600, data_ratio(1,:),'r.','LineWidth',1)
grid on
title('30-minute arc-length')
legend(sprintf('MSIS00/HASDM, mean=%0.2f, rms=%0.2f', nom_mean(1), nom_rms(1)),sprintf('Spire-EDR/HASDM, mean=%0.2f, rms=%0.2f', data_mean(1), data_rms(1)))
set(gca,'FontSize',14)

subplot(3,1,2)
plot(time_rho_all(2,:)/3600, nom_ratio(2,:),'k.','LineWidth',1)
hold on
plot(time_rho_all(2,:)/3600, data_ratio(2,:),'r.','LineWidth',1)
grid on
title('60-minute arc-length')
legend(sprintf('MSIS00/HASDM, mean=%0.2f, rms=%0.2f', nom_mean(2), nom_rms(2)),sprintf('Spire-EDR/HASDM, mean=%0.2f, rms=%0.2f', data_mean(2), data_rms(2)))
set(gca,'FontSize',14)

subplot(3,1,3)
plot(time_rho_all(3,:)/3600, nom_ratio(3,:),'k.','LineWidth',1)
hold on
plot(time_rho_all(3,:)/3600, data_ratio(3,:),'r.','LineWidth',1)
grid on
title('90-minute arc-length')
legend('MSIS00/HASDM','Spire-EDR/HASDM')
set(gca,'FontSize',14)
legend(sprintf('MSIS00/HASDM, mean=%0.2f, rms=%0.2f', nom_mean(3), nom_rms(3)),sprintf('Spire-EDR/HASDM, mean=%0.2f, rms=%0.2f', data_mean(3), data_rms(3)))  
saveas(gcf, fullfile(fig_path, strcat('spire_', year, '_', sat_id, '_density_edr_ratio.png'))) 