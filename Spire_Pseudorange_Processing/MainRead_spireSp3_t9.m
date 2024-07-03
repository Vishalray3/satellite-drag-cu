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
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods/data/ancillary_data')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/PosVelTransformations')
%%%% Files path (rinex, GPS navigation, quaternions & sp3 files):
dir_sp3_q='/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/FM133/leoOrb';

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
yyyy_sat = '2023';
mon_sat = '09';
day_sat = {'05', '06', '07'};
step_size = 60;               % index step at which to store data
%% Read EOP for eci to ecef conversion

% fprintf('Reading Earth Orientation parameters\n');strlen = 37;
% eop = read_finals_data(2017);
% eop.dat = read_tai_utc_dat; % reads only the most current value in the local data file
% eqeterms = 1;


load LeapSeconds.mat
load EOP.mat
load Nutation1980.mat
eqeterms = 1;
EpochUTC = GREGORIANtoJD_vector(str2double(yyyy_sat),str2double(mon_sat),str2double(day_sat{1}));
leap_sec = LeapSeconds(EpochUTC, LeapSecInfo);
%% Find satellite ids
dir_struct = dir(dir_sp3_q);
dir_name_all = {dir_struct.name};
dir_name = [];
for ii = 1:numel(day_sat)
    date_curr(ii,:) = datestr([str2double(yyyy_sat),str2double(mon_sat),str2double(day_sat{ii}),00,00,00],'yyyy-mm-dd');
    dir_name = [dir_name,dir_name_all(contains(dir_name_all, date_curr(ii,:)))];
end

id_common = {'FM133'};
% id_common = strcat('FM',unique(nav_ids));
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
yyyy_sat_num = str2double(yyyy_sat);
mon_sat_num = str2double(mon_sat);
day_sat_num = str2double(day_sat);
data_ = read_sp3_leoOrb(yyyy_sat_num,mon_sat_num,day_sat_num,sat_ID, dir_sp3_q);
% something wrong with the Quaternion data in the version 1.1 file, use
% quaternions from version 1.3. But posistions/velocities seem more
% accurate in vesion 1.1.

%% Remover overlapping data
[data_, time_sec_sp3, time_sec_att, sod_all, pos_error_true, vel_error_true, sod_error_true] = ...
    remove_data_overlap(data_, day_sat_num);

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


for ii = 1:numel(data_.ss)
    %     [reci1,veci1, T_eci2ecef] = ecef2eci  ( data_.sp3_p(:,ii), data_.sp3_v(:,ii),ttt(ii),jdut1(ii),lod(ii),xp(ii),yp(ii),eqeterms,dpsi(ii), deps(ii), nut80);
    %     Xsp3_eci(:,ii) = [reci1;veci1];
    
    [reci1,veci1,rot_ECI2ECEF,dut1,xp,yp] = PosVelConvert(data_.sp3_p(:,ii)'/1000, data_.sp3_v(:,ii)'/1000, jdutc_sec(ii),...
        'ECF2J2K','106terms', leap_sec, NUT1980Info, EOPInfo);
    reci1 = reci1'*1e3;
    veci1 = veci1'*1e3;
    
    Xsp3_eci(:,ii) = [reci1; veci1];
    %     eci_error(:,ii) = Xsp3_eci(:,ii) - Xsp3_eci_cara(:,ii);
    
    
    theta_data(ii) = 0;         % angle in x-y plane
    phi_data(ii) = 0;  % angle with z-axis
    
    %     rot_SBF2ECI(:,:,ii) = rot_ECItoVVLH'*rot_VVLHtoSBF'; % VVLH2ECI*SBF2VVLH
end
rms_pos = rms(pos_error_true, 2);
rms_vel = rms(vel_error_true, 2);

std_pos = std(pos_error_true, 0, 2);
std_vel = std(vel_error_true, 0, 2);

q_sp3 =[];
q_uni_all = [];
yaw_uni = [];
yaw_uni_all = [];
pitch_uni = [];
q_uni =[];
time_sec_all = [];
pitch_uni_all = [];
roll_uni = [];
roll_uni_all = [];
save(strcat('spire_satellite_data_',num2str(sat_ID),'_',yyyy_sat,'_',mon_sat,'_',day_sat{1}),'data_','Xsp3_eci','yyyy_calc','mon_calc','day_calc',...
    'hh_calc','mm_calc','ss_calc','time_sec_sp3','q_sp3','theta_data','phi_data','time_uni','q_uni','time_sec_all', 'sod_error_true', 'pos_error_true',...
    'vel_error_true', 'std_pos', 'std_vel', 'q_uni_all', 'yaw_uni', 'yaw_uni_all', 'pitch_uni', 'pitch_uni_all', 'roll_uni', 'roll_uni_all', 'time_sec_att');
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