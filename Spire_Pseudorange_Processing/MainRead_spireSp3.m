%% Read Spire sp3 and attitude files
% May 30, 2022 :  Vishal Ray
% Colorado Center for Astrodynamics Research
% Copyright 2021 University of Colorado, Boulder
%% Comments
% save data to the OD simulation path
%%
clc
clearvars
restoredefaultpath

linux_os = 0;
if linux_os == 1
    parent_directory = '/home/vira0155';
    dir_sp3_q = '/media/faraday/DATA/thermospheric/spire_data/2022/spire';
else
    parent_directory = '/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder';
    dir_sp3_q='/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Spire 2022 February';
end
addpath(fullfile(parent_directory, 'satellite-drag-cu/Density_inversion_methods/data/ancillary_data'))
addpath(fullfile(parent_directory, 'satellite-drag-cu/PosVelTransformations'))

%% Global constants
global Omega_E c f1 f2 muE
% Omega_E=7.292115*10^-5; %7.2921151467*10^-5; % rad/sec , WGS84 Earth rotation rate
c   = 2.99792458e8;    % GPS acceptd speed of light, m/s
f1  =  1575.42*10^6;   % GPS L1 frequency, in Hz
f2  =  1227.60*10^6;    % GPS L2 frequency, in Hz
muE = 3.986005e14;     % WGS-84 value, m^3/s^2
Omega_E = 7.2921151467*10^-5; %7292115.8553e-11 + 4.3e-15*(GREGORIANtoJD_vector(yyyy,mon,day) - 2451545)/36525;
%% Inputs
% sat_ID = 'FM104';
yyyy_sat = 2022;
mon_sat = 2;
day_sat = [2, 3];
step_size = 1;               % index step at which to store data 
%% Read EOP for eci to ecef conversion

% fprintf('Reading Earth Orientation parameters\n');strlen = 37;
% eop = read_finals_data(2017);
% eop.dat = read_tai_utc_dat; % reads only the most current value in the local data file
% eqeterms = 1;


load LeapSeconds.mat
load EOP.mat
load Nutation1980.mat
eqeterms = 1;
EpochUTC = GREGORIANtoJD_vector(yyyy_sat, mon_sat, day_sat(1));
leap_sec = LeapSeconds(EpochUTC, LeapSecInfo);
%% Find satellite ids
dir_struct = dir(dir_sp3_q);
dir_name_all = {dir_struct.name};
dir_name = [];
for ii = 1:numel(day_sat)
    date_curr(ii,:) = datestr([yyyy_sat, mon_sat, day_sat(ii), 00, 00, 00],'yyyy-mm-dd');
    dir_name = [dir_name,dir_name_all(contains(dir_name_all, date_curr(ii,:)))];
end
telAtt_files = dir_name(contains(dir_name, 'telAtt'));
nav_files = dir_name(contains(dir_name, 'navigation'));
telAtt_ids = extractBetween(telAtt_files, 'FM', '_telAtt');
nav_ids = extractBetween(nav_files, 'FM', '_navigation');
% id_common = strcat('FM',intersect(telAtt_ids, nav_ids));
id_common = {'FM100'};
%% Read data
rms_pos_all = [];
rms_vel_all = [];
for ii = 1:numel(id_common)
    ii
    sat_ID = id_common{ii};
    [rms_pos, rms_vel] = run_datageneration(dir_sp3_q, yyyy_sat,mon_sat,day_sat, sat_ID, leap_sec, NUT1980Info, EOPInfo, step_size);
    if isempty(rms_pos)
        rms_pos = NaN(3,1);
        rms_vel = NaN(3,1);
    end
    rms_pos_all = [rms_pos_all, rms_pos];
    rms_vel_all = [rms_vel_all, rms_vel];
end


function [rms_pos, rms_vel] = run_datageneration(dir_sp3_q, yyyy_sat, mon_sat, day_sat, sat_ID, leap_sec, NUT1980Info, EOPInfo, step_size)
rms_pos = [];
rms_vel = [];

% data_ = read_sp3_leoOrb(yyyy_sat_num,mon_sat_num,day_sat_num,sat_ID,
% dir_sp3_q); % for the new directory structure
data_ = read_sp3_att(yyyy_sat, mon_sat, day_sat, sat_ID, dir_sp3_q);   

% something wrong with the Quaternion data in the version 1.1 file, use
% quaternions from version 1.3. But positions/velocities seem more
% accurate in version 1.1.

%% Remover overlapping data
[data_, time_sec_sp3, time_sec_att, sod_all, pos_error_true, vel_error_true, sod_error_true] = ...
    remove_data_overlap(data_, day_sat);

q_sp3log(1,:) = data_.qx;   q_sp3log(2,:) = data_.qy;      q_sp3log(3,:) = data_.qz;     q_sp3log(4,:) = data_.qw;
%% Complete attitude and orb file
time_sec_all = 86400*GREGORIANtoJD_vector(data_.yyyy_att_all,data_.mon_att_all,data_.dday_att_all) +...
    3600*data_.hh_att_all + 60*data_.mm_att_all + data_.ss_att_all;


%% Used by Vallado
% % Julian Centuries past 1-Jan-2000 12:00 Terrestrial Time
% ttt = (jdutc_sec + eop.dat + 32.184 - 86400*2451545.0)/86400/36525; % sec
% % Interpolate Delta UT1 to TLE time stamps
% dut1 = interp1(86400*eop.fmjd,eop.dut1,jdutc_sec-86400*2400000.5,'linear'); % still in sec
% % Julian day in UT1
% jdut1 = (jdutc_sec+dut1)/86400;
% % Interpolate Earth Orientation parameters (EOP) to TLE time stamps (in UTC)
% lod = 1e-3*interp1(86400*eop.fmjd,eop.rlod,jdutc_sec-86400*2400000.5,'linear'); % now in sec
% xp = (pi/180/3600)*interp1(86400*eop.fmjd,eop.xp,jdutc_sec-86400*2400000.5,'linear'); % now in radians (from arc-sec)
% yp = (pi/180/3600)*interp1(86400*eop.fmjd,eop.yp,jdutc_sec-86400*2400000.5,'linear'); % now in radians (from arc-sec)
% dpsi = (pi/180/3600/1000)*interp1(86400*eop.fmjd,eop.dpsi,jdutc_sec-86400*2400000.5,'linear'); % now in radians (from milli-arc-sec)
% deps = (pi/180/3600/1000)*interp1(86400*eop.fmjd,eop.deps,jdutc_sec-86400*2400000.5,'linear'); % now in radians (from milli-arc-sec)
% % vcorot is for transforming ECEF atmosphere co-rotation v = [0;0;0] and wind = [un;vn;wn] to ECI
% vcorot = [0;0;0];


%% Save the three sets of data: Sp3 eci, calculated eci, log sp3 eci and times corresponding to that

[time_com,ind_com, ind_com2] = intersect(time_sec_att,time_sec_sp3);
if numel(time_com) == numel(time_sec_sp3)
    q_sp3 = q_sp3log(:,ind_com);
else
    q_sp3 = NaN(4, numel(time_sec_sp3));
    q_sp3(:,ind_com2) = q_sp3log(:,ind_com);
end

%% Keep the data at the step size
data_.yyyy = data_.yyyy(1:step_size:end);
data_.mon = data_.mon(1:step_size:end);
data_.dday = data_.dday(1:step_size:end);
data_.hh = data_.hh(1:step_size:end);
data_.mm = data_.mm(1:step_size:end);
data_.ss = data_.ss(1:step_size:end);
data_.sp3_p = data_.sp3_p(:,1:step_size:end);
data_.sp3_v = data_.sp3_v(:,1:step_size:end);
time_sec_sp3 = time_sec_sp3(1:step_size:end);

hh_calc = data_.hh;
mm_calc = data_.mm;
ss_calc = data_.ss;
yyyy_calc = data_.yyyy;
mon_calc = data_.mon;
day_calc = data_.dday;

% Uniformly spaced time vector at every 10 s
time_uni = time_sec_sp3(1):10:time_sec_sp3(end);

%% Convert ECI to ECEF (ITRF) - sp3 file
jdgps_sec = 86400*GREGORIANtoJD_vector(data_.yyyy,data_.mon,data_.dday) + 3600*data_.hh + 60*data_.mm + data_.ss;
jdutc_sec = jdgps_sec - leap_sec + 19;
%% Exit if there's no intersection between leoOrb and leoAtt files
if ~((time_sec_att(1) > time_uni(end)) || (time_uni(1) > time_sec_att(end)))
    
    for ii = 1:numel(data_.ss)
        %     [reci1,veci1, T_eci2ecef] = ecef2eci  ( data_.sp3_p(:,ii), data_.sp3_v(:,ii),ttt(ii),jdut1(ii),lod(ii),xp(ii),yp(ii),eqeterms,dpsi(ii), deps(ii), nut80);
        %     Xsp3_eci(:,ii) = [reci1;veci1];
        
        [reci1,veci1,rot_ECI2ECEF,dut1,xp,yp] = PosVelConvert(data_.sp3_p(:,ii)'/1000, data_.sp3_v(:,ii)'/1000, jdutc_sec(ii),...
            'ECF2J2K','106terms', leap_sec, NUT1980Info, EOPInfo);
        reci1 = reci1'*1e3;
        veci1 = veci1'*1e3;
        
        Xsp3_eci(:,ii) = [reci1; veci1];
        %     eci_error(:,ii) = Xsp3_eci(:,ii) - Xsp3_eci_cara(:,ii);
        
        %% Drag and SRP-related parameters
        rot_ECItoVVLH(3,:) = -Xsp3_eci(1:3,ii)'/norm(Xsp3_eci(1:3,ii));
        yhat               = cross(Xsp3_eci(4:6,ii)/norm(Xsp3_eci(4:6,ii)),Xsp3_eci(1:3,ii)/norm(Xsp3_eci(1:3,ii)));
        rot_ECItoVVLH(2,:) = yhat'/norm(yhat);
        rot_ECItoVVLH(1,:) = cross(rot_ECItoVVLH(2,:),rot_ECItoVVLH(3,:));
        rot_ECItoVVLH(1,:) = rot_ECItoVVLH(1,:)/norm(rot_ECItoVVLH(1,:));
        
        
        %     [reci1,corotwind_eci, T_eci2ecef] = ecef2eci  ( data_.sp3_p(:,ii), [0;0;0],ttt(ii),jdut1(ii),lod(ii),xp(ii),yp(ii),eqeterms,dpsi(ii), deps(ii), nut80);
        [~,corotwind_eci] = PosVelConvert(data_.sp3_p(:,ii)'/1000, [0,0,0],jdutc_sec(ii),'ECF2J2K','106terms', leap_sec, NUT1980Info, EOPInfo);
        corotwind_eci = corotwind_eci'*1e3;
        
        vrel_eci = veci1 - corotwind_eci;
        vrel_vvlh = rot_ECItoVVLH*vrel_eci;
        
        if isnan(q_sp3(1,ii)) && ii > 1
            q_sp3(:,ii) = q_sp3(:,ii-1);
        elseif isnan(q_sp3(:,1))
            ind_notnan = find(isnan(q_sp3(:,1)));
            q_sp3(:,1) = q_sp3(:,ind_notnan(1));
        end
        
        rot_VVLHtoSBF = quat2rotmat(q_sp3(:,ii));
        
        vrel_sbf = rot_VVLHtoSBF*vrel_vvlh;
        
        theta_data(ii) = atan2(vrel_sbf(2),vrel_sbf(1));         % angle in x-y plane
        phi_data(ii) = atan2(vrel_sbf(3),sqrt(vrel_sbf(1)^2+vrel_sbf(2)^2));  % angle with z-axis
        
        %     rot_SBF2ECI(:,:,ii) = rot_ECItoVVLH'*rot_VVLHtoSBF'; % VVLH2ECI*SBF2VVLH
    end
    %% 
    q_temp = quaternion(q_sp3log(4,:), q_sp3log(1,:), q_sp3log(2,:), q_sp3log(3,:));
    q_temp_all = quaternion(data_.qw_all, data_.qx_all, data_.qy_all, data_.qz_all);
    for ii = 1:numel(time_uni)
        if ~isempty(time_sec_att)
            vec_sod = time_sec_att - time_uni(ii);
            ind_sod = find(vec_sod<0, 1, 'last');
            if ind_sod < numel(time_sec_att)
                frac_sec = (time_uni(ii) - time_sec_att(ind_sod))/(time_sec_att(ind_sod+1) - time_sec_att(ind_sod));
                q_interp = slerp(q_temp(:,ind_sod),q_temp(:,ind_sod+1),frac_sec);
                q_mat       = compact(q_interp);
                q_uni(:,ii) = [q_mat(2:4)'; q_mat(1)];
                q_uni(1:3,ii) = q_uni(1:3,ii)./sqrt(1-q_uni(4,ii).^2); % divide by sin(w/2)
                q_uni(1:3,ii) = q_uni(1:3,ii)./sqrt(sum(q_uni(1:3,ii).^2,1)); % normalize by {qx,qy,qz} magnitude
                q_uni(1:3,ii) = q_uni(1:3,ii).*sqrt(1-q_uni(4,ii).^2); % multiply by sin(w/2)
            elseif isempty(ind_sod)
                q_uni(:,ii) = q_sp3log(:,1);
            else
                q_uni(:,ii) = q_uni(:,ii-1);
            end
        else
            q_uni = [];
            yaw_uni = [];
            pitch_uni = [];
            roll_uni = [];
        end
        
        if ~isempty(time_sec_all)
            vec_sod = time_sec_all - time_uni(ii);
            ind_sod = find(vec_sod<0, 1, 'last');
            if ind_sod < numel(time_sec_all)
                frac_sec = (time_uni(ii) - time_sec_all(ind_sod))/(time_sec_all(ind_sod+1) - time_sec_all(ind_sod));
                q_interp = slerp(q_temp_all(:,ind_sod),q_temp_all(:,ind_sod+1),frac_sec);
                q_mat       = compact(q_interp);
                q_uni_all(:,ii) = [q_mat(2:4)'; q_mat(1)];
                q_uni_all(1:3,ii) = q_uni_all(1:3,ii)./sqrt(1-q_uni_all(4,ii).^2); % divide by sin(w/2)
                q_uni_all(1:3,ii) = q_uni_all(1:3,ii)./sqrt(sum(q_uni_all(1:3,ii).^2,1)); % normalize by {qx,qy,qz} magnitude
                q_uni_all(1:3,ii) = q_uni_all(1:3,ii).*sqrt(1-q_uni_all(4,ii).^2); % multiply by sin(w/2)
            elseif isempty(ind_sod)
                q_uni_all(:,ii) = [data_.qx_all(1);data_.qy_all(1);data_.qz_all(1);data_.qw_all(1)];
            else
                q_uni_all(:,ii) = q_uni_all(:,ii-1);
            end
        else
            q_uni_all = [];
        end
    end
    
    if ~isempty(time_sec_att)  
        [yaw_uni, pitch_uni, roll_uni] = quat2angle([q_uni(4,:)',q_uni(1,:)', q_uni(2,:)', q_uni(3,:)']);
    else
        yaw_uni = [];
        pitch_uni = [];
        roll_uni = [];
    end
    if ~isempty(time_sec_all)
        [yaw_uni_all, pitch_uni_all, roll_uni_all] = quat2angle([q_uni_all(4,:)',q_uni_all(1,:)', q_uni_all(2,:)', q_uni_all(3,:)']);
    else
        yaw_uni_all = [];
        pitch_uni_all = [];
        roll_uni_all = [];        
    end
    
%     [yaw_uni_data, pitch_uni_data, roll_uni_data] = quat2angle([q_sp3log(4,:)',q_sp3log(1,:)', q_sp3log(2,:)', q_sp3log(3,:)']);
%     [yaw_uni_all_data, pitch_uni_all_data, roll_uni_all_data] = quat2angle([data_.qw_all',data_.qx_all',data_.qy_all', data_.qz_all']);
    
    
    %%
    %%
    
%     figure()
%     subplot(2,1,1)
%     plot(sod_error_true/3600, abs(pos_error_true(1,:))*1e2,'.')
%     hold on
%     plot(sod_error_true/3600, abs(pos_error_true(2,:))*1e2,'.')
%     plot(sod_error_true/3600, abs(pos_error_true(3,:))*1e2,'.')
%     ylabel('Position (cm)')
%     legend('X','Y','Z')
%     title('POD overlap errors in ECEF (Spire)')
%     set(gca,'FontSize',16)
%     grid on
%     subplot(2,1,2)
%     plot(sod_error_true/3600, abs(vel_error_true(1,:))*1e3,'.')
%     hold on
%     plot(sod_error_true/3600, abs(vel_error_true(2,:))*1e3,'.')
%     plot(sod_error_true/3600, abs(vel_error_true(3,:))*1e3,'.')
%     ylabel('Velocity (mm/s)')
%     set(gca,'FontSize',16)
%     grid on
%     xlabel('Hours')
%     
%     sod_uni = time_uni;
%     sod_att = time_sec_att;
%     sod_att_all = time_sec_all;
%     figure(3)
%     subplot(3,1,1)
%     plot(sod_uni/3600, yaw_uni*180/pi,'.')
%     hold on
%     plot(sod_uni/3600, yaw_uni_all*180/pi,'.')
%     plot(sod_att/3600, yaw_uni_data*180/pi,'.')
%     plot(sod_att_all/3600, yaw_uni_all_data*180/pi,'.')
%     ylabel('Yaw (deg)')
%     legend('LeoAtt (high-cadence, data-gaps)', 'TelAtt (low-cadence, continuous)', 'LeoAtt data', 'TelAtt data')
%     title('Attitude data')
%     set(gca,'FontSize',16)
%     grid on
%     subplot(3,1,2)
%     plot(sod_uni/3600, pitch_uni*180/pi,'.')
%     hold on
%     plot(sod_uni/3600, pitch_uni_all*180/pi,'.')
%     plot(sod_att/3600, pitch_uni_data*180/pi,'.')
%     plot(sod_att_all/3600, pitch_uni_all_data*180/pi,'.')
%     ylabel('Pitch (deg)')
%     title('Attitude data')
%     set(gca,'FontSize',16)
%     grid on
%     subplot(3,1,3)
%     plot(sod_uni/3600, roll_uni*180/pi,'.')
%     hold on
%     plot(sod_uni/3600, roll_uni_all*180/pi,'.')
%     plot(sod_att/3600, roll_uni_data*180/pi,'.')
%     plot(sod_att_all/3600, roll_uni_all_data*180/pi,'.')
%     ylabel('Roll (deg)')
%     title('Attitude data')
%     set(gca,'FontSize',16)
%     grid on
%     xlabel('Hours')
    
    rms_pos = rms(pos_error_true, 2);
    rms_vel = rms(vel_error_true, 2);
    
    std_pos = std(pos_error_true, 0, 2);
    std_vel = std(vel_error_true, 0, 2);
    
    
    save(strcat('spire_satellite_data_',sat_ID,'_',num2str(yyyy_sat),'_',sprintf('%02d', mon_sat),'_', sprintf('%02d', day_sat(1))),'data_','Xsp3_eci','yyyy_calc','mon_calc','day_calc',...
        'hh_calc','mm_calc','ss_calc','time_sec_sp3','q_sp3','theta_data','phi_data','time_uni','q_uni','time_sec_all', 'sod_error_true', 'pos_error_true',...
        'vel_error_true', 'std_pos', 'std_vel', 'q_uni_all', 'yaw_uni', 'yaw_uni_all', 'pitch_uni', 'pitch_uni_all', 'roll_uni', 'roll_uni_all', 'time_sec_att');
end
end
%% Plots overlapping data

% pos_diff = pos_all(:,ind_over) - pos_all(:,ind_data_common);
% vel_diff = vel_all(:,ind_over) - vel_all(:,ind_data_common);
% figure(1)
% subplot(2,1,1)
% plot(sod_all(ind_over)/3600, abs(pos_diff(1,:))*1e2,'.')
% hold on
% plot(sod_all(ind_over)/3600, abs(pos_diff(2,:))*1e2,'.')
% plot(sod_all(ind_over)/3600, abs(pos_diff(3,:))*1e2,'.')
% ylabel('Position (cm)')
% legend('X','Y','Z')
% title('POD overlap errors in ECEF')
% set(gca,'FontSize',16)
% grid on
% subplot(2,1,2)
% plot(sod_all(ind_over)/3600, abs(vel_diff(1,:))*1e3,'.')
% hold on
% plot(sod_all(ind_over)/3600, abs(vel_diff(2,:))*1e3,'.')
% plot(sod_all(ind_over)/3600, abs(vel_diff(3,:))*1e3,'.')
% ylabel('Velocity (mm/s)')
% set(gca,'FontSize',16)
% grid on
% xlabel('Hours')