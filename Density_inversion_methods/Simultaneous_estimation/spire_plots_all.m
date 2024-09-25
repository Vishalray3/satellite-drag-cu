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


sat_ids_skip = {};
year_data = 2022;
%%
dir_struct = dir(output_dir);
dir_name = {dir_struct.name};
data_pattern = "results_cdpanel_spire_satellite_data_FM" + digitsPattern(3) + '_' + digitsPattern(4) + '_' + digitsPattern(2) + '_' + digitsPattern(2) + '.mat';

file_names_all = dir_name(matches(dir_name, data_pattern));

sat_ID_mat = unique(extractBetween(file_names_all, strcat('data', data_pod, '_'),  strcat('_', sprintf('%d', year_data)))); %  %
dataset_sat_mat = extractBefore(file_names_all, '_satellite');
sat_ID_mat(ismember(sat_ID_mat, sat_ids_skip)) = [];

%% Run the loop

for sat_id = sat_ID_mat
    rho_data_all = [];
    rho_nom_all = [];
    rho_hasdm_all = [];
    time_rho_all = [];
    file_names_idx = find(contains(file_names_all, sat_id));
    file_names = file_names_all(file_names_idx);
    for ii=1:numel(file_names)
        load(fullfile(output_dir, file_names{ii}))
        rho_hasdm_all = [rho_hasdm_all, rho_hasdm_eff];
        rho_nom_all = [rho_nom_all, rho_nom_eff];
        rho_data_all = [rho_data_all, rho_eff];
        time_cal = datetime((jd_init + time_rho)/86400,'convertfrom','juliandate','Format','yyy-MM-dd-hh-mm-ss');
        time_rho_all = [time_rho_all, time_cal];
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
    plot(time_rho_all(1,:),rho_hasdm_all(1,:),'k.','LineWidth',1)
    hold on
    plot(time_rho_all(1,:),rho_nom_all(1,:),'r.','LineWidth',1)
    plot(time_rho_all(1,:),rho_data_all(1,:),'g.','LineWidth',1)
    grid on
    title('30-minute arc-length')
    legend('HASDM','MSIS00','Spire EDR')
    set(gca,'FontSize',14)

    subplot(3,1,2)
    plot(time_rho_all(2,:),rho_hasdm_all(2,:),'k.','LineWidth',1)
    hold on
    plot(time_rho_all(2,:),rho_nom_all(2,:),'r.','LineWidth',1)
    plot(time_rho_all(2,:),rho_data_all(2,:),'g.','LineWidth',1)
    grid on
    ylabel('Density ($kg/m^3$)','Interpreter','latex')
    title('60-minute arc-length')
    set(gca,'FontSize',14)

    subplot(3,1,3)
    plot(time_rho_all(3,:),rho_hasdm_all(3,:),'k.','LineWidth',1)
    hold on
    plot(time_rho_all(3,:),rho_nom_all(3,:),'r.','LineWidth',1)
    plot(time_rho_all(3,:),rho_data_all(3,:),'g.','LineWidth',1)
    grid on
    ylabel('Density ($kg/m^3$)','Interpreter','latex')
    title('90-minute arc-length')
    set(gca,'FontSize',14)
    saveas(gcf, fullfile(fig_path, strcat('spire_', sprintf('%d', year_data), '_', sat_id, '_density_edr.png')))    

    figure(2)
    subplot(3,1,1)
    plot(time_rho_all(1,:), nom_ratio(1,:),'k.','LineWidth',1)
    hold on
    plot(time_rho_all(1,:), data_ratio(1,:),'r.','LineWidth',1)
    grid on
    title('30-minute arc-length')
    legend(sprintf('MSIS00/HASDM, mean=%0.2f, rms=%0.2f', nom_mean(1), nom_rms(1)),sprintf('Spire-EDR/HASDM, mean=%0.2f, rms=%0.2f', data_mean(1), data_rms(1)))
    set(gca,'FontSize',14)

    subplot(3,1,2)
    plot(time_rho_all(2,:), nom_ratio(2,:),'k.','LineWidth',1)
    hold on
    plot(time_rho_all(2,:), data_ratio(2,:),'r.','LineWidth',1)
    grid on
    title('60-minute arc-length')
    legend(sprintf('MSIS00/HASDM, mean=%0.2f, rms=%0.2f', nom_mean(2), nom_rms(2)),sprintf('Spire-EDR/HASDM, mean=%0.2f, rms=%0.2f', data_mean(2), data_rms(2)))
    set(gca,'FontSize',14)

    subplot(3,1,3)
    plot(time_rho_all(3,:), nom_ratio(3,:),'k.','LineWidth',1)
    hold on
    plot(time_rho_all(3,:), data_ratio(3,:),'r.','LineWidth',1)
    grid on
    title('90-minute arc-length')
    legend('MSIS00/HASDM','Spire-EDR/HASDM')
    set(gca,'FontSize',14)
    legend(sprintf('MSIS00/HASDM, mean=%0.2f, rms=%0.2f', nom_mean(3), nom_rms(3)),sprintf('Spire-EDR/HASDM, mean=%0.2f, rms=%0.2f', data_mean(3), data_rms(3)))  
    saveas(gcf, fullfile(fig_path, strcat('spire_', sprintf('%d', year_data), '_', sat_id, '_density_edr_ratio.png')))
end