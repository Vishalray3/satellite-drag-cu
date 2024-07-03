%% compare ucar and spire pods

clc
clearvars

load ucar_02_05_2022

rms_pos_all_ucar = rms_pos_all;
rms_vel_all_ucar = rms_vel_all;

load spire_02_05_2022

rms_pos_all_spire = rms_pos_all;
rms_vel_all_spire = rms_vel_all;

id_all = erase(id_all,'FM');
[union_ids, ind_ucar, ind_spire] = union(ucar_id_all, id_all);

rms_ucar = NaN(2, numel(union_ids));
rms_spire = NaN(2, numel(union_ids));
for ii = 1:numel(union_ids)
    ind_ucar = strcmp(ucar_id_all,union_ids(ii));
    if any(ind_ucar)
        rms_ucar(:, ii) = [vecnorm(rms_pos_all_ucar(:, ind_ucar), 2,1); vecnorm(rms_vel_all_ucar(:, ind_ucar),2,1)];
    end
    
    ind_spire = strcmp(id_all, union_ids(ii));
    if any(ind_spire)
        rms_spire(:, ii) = [vecnorm(rms_pos_all_spire(:, ind_spire),2,1); vecnorm(rms_vel_all_spire(:, ind_spire),2,1)];
    end
end

%% Plot bar charts
rms_pos = [rms_ucar(1,:)', rms_spire(1,:)']; 
rms_vel = [rms_ucar(2,:)', rms_spire(2,:)']; 

barx = categorical(union_ids);
subplot(2,1,1)
bar(barx,rms_pos*1e2, 1)
title('RMS of orbit overlaps for Spire and UCAR PODs')
ylabel('RMS position (cm)')
set(gca,'FontSize',14)
grid on
legend('UCAR', 'Spire')
subplot(2,1,2)
bar(barx,rms_vel*1e3)
ylabel('RMS velocity (mm/s)')
xlabel('Satellite IDs')
set(gca,'FontSize',14)
grid on

%%    
common_ids = intersect(ucar_id_all, id_all);
rms_pos = [];
rms_vel = [];
for ii = 1:numel(common_ids)
load(strcat('spire_satellite_data_FM', common_ids{ii}, '_2022_02_05'))
data_spire = data_;
time_spire = time_sec_sp3;

load(strcat('spire_satellite_data_ucar_FM', common_ids{ii}, '_2022_02_05'))
data_ucar = data_;
time_ucar = time_sec_sp3;

[time_spire, ind_spire, ind_ucar] = intersect(time_spire, time_ucar);

rms_pos(ii) = norm(rms(data_spire.sp3_p(:, ind_spire) - data_ucar.sp3_p(:, ind_ucar), 2));
rms_vel(ii) = norm(rms(data_spire.sp3_v(:, ind_spire) - data_ucar.sp3_v(:, ind_ucar), 2));
end


barx = categorical(common_ids);
subplot(2,1,1)
bar(barx,rms_pos*1e2, 1)
title('RMS of differences between Spire and UCAR PODs')
ylabel('RMS position (cm)')
set(gca,'FontSize',14)
grid on

subplot(2,1,2)
bar(barx,rms_vel*1e3)
ylabel('RMS velocity (mm/s)')
xlabel('Satellite IDs')
set(gca,'FontSize',14)
grid on

% sod_error_true = time_spire - time_sec_sp3(1);
% figure(2)
% subplot(2,1,1)
% plot(sod_error_true/3600, abs(pos_error_true(1,:))*1e2,'.')
% hold on
% plot(sod_error_true/3600, abs(pos_error_true(2,:))*1e2,'.')
% plot(sod_error_true/3600, abs(pos_error_true(3,:))*1e2,'.')
% ylabel('Position (cm)')
% legend('X','Y','Z')
% title('Differences between Spire and UCAR')
% set(gca,'FontSize',16)
% grid on
% subplot(2,1,2)
% plot(sod_error_true/3600, abs(vel_error_true(1,:))*1e3,'.')
% hold on
% plot(sod_error_true/3600, abs(vel_error_true(2,:))*1e3,'.')
% plot(sod_error_true/3600, abs(vel_error_true(3,:))*1e3,'.')
% ylabel('Velocity (mm/s)')
% set(gca,'FontSize',16)
% grid on
% xlabel('Hours')