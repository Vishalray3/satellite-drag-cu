%% Positions and velocities from raw pseudorange measurements
% Kinematic solution from raw pseudorange measurements for the SPIRE
% satellites
% Aug 05, 2021 :  Vishal Ray
% Colorado Center for Astrodynamics Research
% Copyright 2021 University of Colorado, Boulder
clc
clearvars
restoredefaultpath
%%%% Files path (rinex, GPS navigation, quaternions & sp3 files):
dir_rnx = '/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Density inversion backup/Spire_Pseudorange_Processing/DataFiles/spire_20180923_20181209_antPOD_83-85/';
dir_sp3_q='/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Density inversion backup/Spire_Pseudorange_Processing/DataFiles/spire_20180923_20181209_OrbAtt_83-85/';
dir_GPSn='/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Density inversion backup/Spire_Pseudorange_Processing/DataFiles/GPS_nav/';

%% Global constants
global Omega_E c f1 f2 muE eop nut80 eqeterms
% Omega_E=7.292115*10^-5; %7.2921151467*10^-5; % rad/sec , WGS84 Earth rotation rate
c   =  2.99792458e8;    % GPS acceptd speed of light, m/s
f1  =  1575.42*10^6;    % GPS L1 frequency, in Hz
f2  =  1227.60*10^6;    % GPS L2 frequency, in Hz
muE =  3.986005e14;     % WGS-84 value, m^3/s^2
%% Inputs
load('spire_leoAtt_083', 'reci', 'veci', 'doy', 'hh', 'mm', 'ss', 'sod', 'yyyy', 'theta','mon','day','Cd','Aref', 'rot_SBF2ECI', 'rot_ECI2ECEF','q')
hh_pre = hh;
mm_pre = mm;
ss_pre = ss;
yyyy_pre = yyyy;
mon_pre = mon;
day_pre = day;
rot_SBF2ECI_pre = rot_SBF2ECI;
theta_pre = theta;
q_pre = q;
Cd_pre = Cd;
Aref_pre = Aref;

sat_ID = 85;
pos_iter = 10;
vel_iter = 10;
flag_w = 0;

%% Read EOP for eci to ecef conversion
load nut80.dat
fprintf('Reading Earth Orientation parameters\n');strlen = 37;
eop = read_finals_data;
eop.dat = read_tai_utc_dat; % reads only the most current value in the local data file
eqeterms = 1;
%% Read GPS data and compute all required states
yyyy = 2018;
mon = 11;
day = 7;
Omega_E = 7.2921151467*10^-5; %7292115.8553e-11 + 4.3e-15*(GREGORIANtoJD_vector(yyyy,mon,day) - 2451545)/36525;
[gps_ephem] = read_gpspos_sp3('nga20263');
% [data_] = Spire_Pseudorange_Processing_sp3(yyyy,mon,day,sat_ID,pos_iter,vel_iter,flag_w, gps_ephem, dir_rnx, dir_sp3_q, dir_GPSn);
[data_] = Spire_Pseudorange_Processing_v2(yyyy,mon,day,sat_ID,pos_iter,vel_iter,flag_w, dir_rnx, dir_sp3_q, dir_GPSn);
% something wrong with the Quaternion data in the version 1.1 file, use
% quaternions from version 1.3. But posistions/velocities seem more
% accurate in vesion 1.1.


%% Convert ECI to ECEF (ITRF)
jdgps_sec = 86400*GREGORIANtoJD_vector(yyyy,mon,day) + 3600*data_.hh + 60*data_.mm + data_.ss;
jdutc_sec = jdgps_sec - eop.dat + 19;
% Julian Centuries past 1-Jan-2000 12:00 Terrestrial Time
ttt = (jdutc_sec + eop.dat + 32.184 - 86400*2451545.0)/86400/36525; % sec
% Interpolate Delta UT1 to TLE time stamps
dut1 = interp1(86400*eop.fmjd,eop.dut1,jdutc_sec-86400*2400000.5,'linear'); % still in sec
% Julian day in UT1
jdut1 = (jdutc_sec+dut1)/86400;
% Interpolate Earth Orientation parameters (EOP) to TLE time stamps (in UTC)
lod = 1e-3*interp1(86400*eop.fmjd,eop.rlod,jdutc_sec-86400*2400000.5,'linear'); % now in sec
xp = (pi/180/3600)*interp1(86400*eop.fmjd,eop.xp,jdutc_sec-86400*2400000.5,'linear'); % now in radians (from arc-sec)
yp = (pi/180/3600)*interp1(86400*eop.fmjd,eop.yp,jdutc_sec-86400*2400000.5,'linear'); % now in radians (from arc-sec)
dpsi = (pi/180/3600/1000)*interp1(86400*eop.fmjd,eop.dpsi,jdutc_sec-86400*2400000.5,'linear'); % now in radians (from milli-arc-sec)
deps = (pi/180/3600/1000)*interp1(86400*eop.fmjd,eop.deps,jdutc_sec-86400*2400000.5,'linear'); % now in radians (from milli-arc-sec)
% vcorot is for transforming ECEF atmosphere co-rotation v = [0;0;0] and wind = [un;vn;wn] to ECI
vcorot = [0;0;0];
%% Find common data with the sp3 log files
ind_mon = find(mon_pre == mon);
ind_day = find(day_pre == day);
ind_com = intersect(ind_mon, ind_day);
X_sp3log = [reci(:,ind_com);veci(:,ind_com)]*1e3;
jdgps_sec_sp3     = 86400*GREGORIANtoJD_vector(yyyy,mon,day) + 3600*hh_pre(ind_com) + 60*mm_pre(ind_com) + ss_pre(ind_com);
sod_sec_sp3log = jdgps_sec_sp3 - 86400*GREGORIANtoJD_vector(yyyy,mon,day);
sod_sec_sp3 = jdgps_sec - 86400*GREGORIANtoJD_vector(yyyy,mon,day);

% interpolation is not a good idea, the data is very sparse.
% Xsp3_com = interp1(sod_sec, Xsp3_eci', sod_sec_sp3, 'pchip');
% err_sp3 = Xsp3_com'-Xsp3_proc;

% Remove repeating data for interpolation
[sod_sec_sp3log,indm] = unique(sod_sec_sp3log);
X_sp3log = X_sp3log(indm);

%% Save the three sets of data: Sp3 eci, calculated eci, log sp3 eci and times corresponding to that
hh_calc = data_.hh;
mm_calc = data_.mm;
ss_calc = data_.ss;
yyyy_calc = data_.yyyy;
mon_calc = data_.mon;
day_calc = data_.dday;

hh_sp3log   = hh_pre(ind_com);
hh_sp3log = hh_sp3log(indm);
mm_sp3log   = mm_pre(ind_com);
mm_sp3log = mm_sp3log(indm);
ss_sp3log   = ss_pre(ind_com);
ss_sp3log = ss_sp3log(indm);
yyyy_sp3log = yyyy_pre(ind_com);
yyyy_sp3log = yyyy_sp3log(indm);
mon_sp3log  = mon_pre(ind_com);
mon_sp3log = mon_sp3log(indm);
day_sp3log  = day_pre(ind_com);
day_sp3log = day_sp3log(indm);
doy_vec = doy(ind_com);
doy_vec = doy_vec(indm);
theta_sp3log = theta_pre(ind_com);
theta_sp3log = theta_sp3log(indm);
rot_SBF2ECI_sp3log = rot_SBF2ECI_pre(:,:,ind_com);
rot_SBF2ECI_sp3log = rot_SBF2ECI_sp3log(:,:,indm);
q_sp3log    = q_pre(:, ind_com);
q_sp3log = q_sp3log(:,indm);
Cd_eric    = Cd_pre(:, ind_com);
Cd_eric = Cd_eric(:,indm);
Aref_eric    = Aref_pre(:, ind_com);
Aref_eric = Aref_eric(:,indm);

% interpolation is not a good idea, the data is very sparse.
% not used anywhere
theta_interplog = interp1(sod_sec_sp3log, theta_sp3log, sod_sec_sp3, 'linear');
Cd_interplog = interp1(sod_sec_sp3log, Cd_eric, sod_sec_sp3, 'linear');
Aref_interplog = interp1(sod_sec_sp3log, Aref_eric, sod_sec_sp3, 'linear');

% Uniformly spaced time vector
sod_uni = sod_sec_sp3(1):1:sod_sec_sp3log(end);
for ii = 1:numel(data_.ss)
    [reci1,veci1, T_eci2ecef] = ecef2eci  ( data_.recef(:,ii), data_.vecef(:,ii),ttt(ii),jdut1(ii),lod(ii),xp(ii),yp(ii),eqeterms,dpsi(ii), deps(ii),nut80 );
    X_eci(:,ii) = [reci1;veci1];
    [reci1,veci1, T_eci2ecef] = ecef2eci  ( data_.sp3_p(:,ii), data_.sp3_v(:,ii),ttt(ii),jdut1(ii),lod(ii),xp(ii),yp(ii),eqeterms,dpsi(ii), deps(ii), nut80  );
    Xsp3_eci(:,ii) = [reci1;veci1];
    
%     %% Drag and SRP-related parameters
%     rot_ECItoVVLH(3,:) = -Xsp3_eci(1:3,ii)'/norm(Xsp3_eci(1:3,ii));
%     yhat               = cross(Xsp3_eci(4:6,ii)/norm(Xsp3_eci(4:6,ii)),Xsp3_eci(1:3,ii)/norm(Xsp3_eci(1:3,ii)));
%     rot_ECItoVVLH(2,:) = yhat'/norm(yhat);
%     rot_ECItoVVLH(1,:) = cross(rot_ECItoVVLH(3,:),rot_ECItoVVLH(2,:));
%     rot_ECItoVVLH(1,:) = rot_ECItoVVLH(1,:)/norm(rot_ECItoVVLH(1,:));
%     
%     [reci1,corotwind_eci, T_eci2ecef] = ecef2eci  ( data_.sp3_p(:,ii), [0;0;0],ttt(ii),jdut1(ii),lod(ii),xp(ii),yp(ii),eqeterms,dpsi(ii), deps(ii)  );
%     vrel_eci = veci1 - corotwind_eci;
%     vrel_vvlh = rot_ECItoVVLH*vrel_eci;
%     
%    
%     vec_sod = sod_sec_sp3log - sod_sec_sp3(ii);
%     ind_sod = find(vec_sod<0, 1, 'last');
%     frac_sec = (sod_sec_sp3(ii) - sod_sec_sp3log(ind_sod))/(sod_sec_sp3log(ind_sod+1) - sod_sec_sp3log(ind_sod));
% 
%     rot_VVLHtoSBF = quat2rotmat(q(:,ii));
%     
%     vrel_sbf = rot_VVLHtoSBF*vrel_vvlh;
%     
%     theta_data(ii) = atan2(vrel_sbf(2),vrel_sbf(1));         % angle in x-y plane
%     phi_sp3(ii) = atan2(vrel_sbf(3),sqrt(vrel_sbf(1)^2+vrel_sbf(2)^2));  % angle with z-axis
%     
%     rot_SBF2ECI_sp3(:,:,ii) = (rot_ECItoVVLH*rot_VVLHtoSBF)';
end
%%
q_temp = quaternion(q_sp3log(4,:), q_sp3log(1,:), q_sp3log(2,:), q_sp3log(3,:));
for ii = 1:numel(sod_uni)
    vec_sod = sod_sec_sp3log - sod_uni(ii);
    ind_sod = find(vec_sod<0, 1, 'last');
    frac_sec = (sod_uni(ii) - sod_sec_sp3log(ind_sod))/(sod_sec_sp3log(ind_sod+1) - sod_sec_sp3log(ind_sod));
    q_interp = slerp(q_temp(:,ind_sod),q_temp(:,ind_sod+1),frac_sec);
    q_mat       = compact(q_interp);
    q_sp3(:,ii) = [q_mat(2:4)'; q_mat(1)];
    q_sp3(1:3,ii) = q_sp3(1:3,ii)./sqrt(1-q_sp3(4,ii).^2); % divide by sin(w/2)
    q_sp3(1:3,ii) = q_sp3(1:3,ii)./sqrt(sum(q_sp3(1:3,ii).^2,1)); % normalize by {qx,qy,qz} magnitude
    q_sp3(1:3,ii) = q_sp3(1:3,ii).*sqrt(1-q_sp3(4,ii).^2); % multiply by sin(w/2)   
end
%% Frame conversion using SPICE: DOES NOT WORK WELL - USES ITRF93.
% jdgps_epoch   = 86400*GREGORIANtoJD_vector(2000,1,1) + 3600*12 + 60*0 + 0;
% hh_calc = data_.hh;
% mm_calc = data_.mm;
% ss_calc = data_.ss;
% jdgps_sec     = 86400*GREGORIANtoJD_vector(yyyy,mon,day) + 3600*hh_calc + 60*mm_calc + ss_calc;
% tai_time_gps  = jdgps_sec-jdgps_epoch + 19;
% et_time_gps   = cspice_unitim(tai_time_gps, 'TAI', 'TDB');
% time_prop     = et_time_gps;
% time_prop_utc = time_prop - time_prop(1);
%
% for ii = 1:numel(data_.ss)
%     rot_ECI2ECEF = cspice_pxform( 'J2000', 'ITRF93', time_prop(ii) );
%     X_ecef = [data_.recef(:,ii); data_.vecef(:,ii)];
%     Xsp3_ecef = [data_.sp3_p(:,ii); data_.sp3_v(:,ii)];
%     dXsp3_ecef = [data_.dsp3_p(:,ii); data_.dsp3_v(:,ii)];
%     X_eci(1:3,ii) = rot_ECI2ECEF'*X_ecef(1:3);
%     X_eci(4:6,ii) = rot_ECI2ECEF'*X_ecef(4:6);
%     Xsp3_eci(:,ii) = rot_ECI2ECEF'*Xsp3_ecef;
%     dXsp3_eci(:,ii) = rot_ECI2ECEF'*dXsp3_ecef;
% end
%


% dsp3 = vecnorm(data_.dsp3_p(1:3,:),2,1);
% B = lowpass(dsp3, 1e-3); 
% plot(dsp3)
% hold on
% plot(B,'x')
%%


save(strcat('spire_sat',num2str(sat_ID),'_',num2str(yyyy),'_',num2str(mon),'_',num2str(day)),'data_','X_eci','Xsp3_eci','yyyy_calc','mon_calc','day_calc',...
    'hh_calc','mm_calc','ss_calc','X_sp3log','yyyy_sp3log', 'mon_sp3log','day_sp3log','hh_sp3log','mm_sp3log','ss_sp3log','sod_sec_sp3','sod_sec_sp3log',...
    'q_sp3','q_sp3log','theta_sp3log','doy_vec','sod_uni','theta_interplog', 'Cd_eric', 'Aref_eric', 'Cd_interplog','Aref_interplog');
