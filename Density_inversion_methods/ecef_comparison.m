%% Comparing MICE and Vallado frame conversions.
restoredefaultpath
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/DataTools/mice/mice/src/mice/')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/DataTools/mice/mice/lib/' )
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/public/quaternions_example_Spire_shaylah')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/public/quaternions_example_Spire_shaylah/FundamentalsOfAstrodynamicsAndApplications/matlab')
clc
clear all
% load spire_leoAtt_083
fprintf('Reading Earth Orientation parameters\n');strlen = 37;
eop = read_finals_data;
eop.dat = read_tai_utc_dat; % reads only the most current value in the local data file
eqeterms = 1;
%% Furnishing spice kernels to use the spice functions
cspice_furnsh('de430.bsp')                       %% Planetary ephemerides
cspice_furnsh('naif0012.tls')                    %% time leap seconds
cspice_furnsh('earth_assoc_itrf93.tf')           %% Earth ITRF frame
cspice_furnsh('pck00010.tpc')
cspice_furnsh('gm_de431.tpc')                    %% GM values
cspice_furnsh('earth_000101_181111_180820.bpc')

%% Gravity
load GravCoeff_EGM2008
mu_e    = cspice_bodvrd('EARTH', 'GM',3)*1e9;           %% mu of Earth in m3/s2
mu_s    = cspice_bodvrd('SUN', 'GM',3)*1e9;             %% mu of Sun in m3/s2
mu_m    = cspice_bodvrd('MOON', 'GM',3)*1e9;            %% mu of Moon in m3/s2
Rad_e   = cspice_bodvrd('EARTH', 'RADII',3);            %% radius of Earth from the geophysical kernel (need for shadow function of SRP)
Rad_s   = cspice_bodvrd('SUN', 'RADII',3);              %% radius of Earth from the geophysical kernel (need for shadow function of SRP)
Re      = Rad_e(1)*1e3;
Rs      = Rad_s(1)*1e3;
Psun    = 4.56e-6;                                      %% Radiation pressure at 1 au (Doornbos PhD) N/m2
AU      = 149597870660;                                 %% AU value m
omega_e = 7.2921150e-5;  
deg_grav = 120 ;                                         %% degree of spherical harmonics
ord_grav = 120; 
myHighGrav       = HighGrav1;
myHighGrav.mu_e  = mu_e;
myHighGrav.deg   = deg_grav;
myHighGrav.ord   = ord_grav;
myHighGrav.Re    = Re;
myHighGrav.Cbar  = Cbar(1:deg_grav+1, 1:deg_grav+1);
myHighGrav.Sbar  = Sbar(1:deg_grav+1, 1:deg_grav+1);
%% Vallado
load spire_leoAtt_083
jdgps_sec = 86400*GREGORIANtoJD_vector(yyyy,mon,day) + 3600*hh + 60*mm + ss;
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
for i = 1:length(jdutc_sec)
    % Calculate recef (km) and vecef (km/s) from reci (km) and veci (km/s)
    [recef_val(:,i),vecef_val(:,i),~, TEI] = eci2ecef( reci(:,i), veci(:,i), [0;0;0], ttt(i), jdut1(i), lod(i), ...
        xp(i), yp(i), eqeterms, dpsi(i), deps(i) );
    myHighGrav.TEI = TEI;
    accel(:,i) = myHighGrav.compAccelPar(reci(:,i)*1e3, 0);
    
end

%% MICE
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/DataTools/mice/mice/src/mice/')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/DataTools/mice/mice/lib/' )
% Furnishing spice kernels to use the spice functions
cspice_furnsh('de430.bsp')                       %% Planetary ephemerides
cspice_furnsh('naif0012.tls')                    %% time leap seconds
cspice_furnsh('earth_assoc_itrf93.tf')           %% Earth ITRF frame
cspice_furnsh('pck00010.tpc')
cspice_furnsh('gm_de431.tpc')                    %% GM values
cspice_furnsh('earth_000101_181111_180820.bpc')
myHighGrav       = HighGrav;
myHighGrav.mu_e  = mu_e;
myHighGrav.deg   = deg_grav;
myHighGrav.ord   = ord_grav;
myHighGrav.Re    = Re;
myHighGrav.Cbar  = Cbar(1:deg_grav+1, 1:deg_grav+1);
myHighGrav.Sbar  = Sbar(1:deg_grav+1, 1:deg_grav+1);
% Time conversions  (ET epoch - Jan 1, 2000 12:00:00 TDB)
jdgps_epoch = 86400*GREGORIANtoJD_vector(2000,1,1) + 3600*12 + 60*0 + 0;

tai_time_gps  = jdgps_sec-jdgps_epoch + 19;
et_time_gps   = cspice_unitim(tai_time_gps, 'TAI', 'TDB');
time_prop     = et_time_gps;

for ii =1:numel(time_prop)
TEI = cspice_sxform( 'J2000', 'ITRF93', time_prop(ii) ); 
X_ecef(:,ii) = TEI*[reci(:,ii);veci(:,ii)];
myHighGrav.et = time_prop(ii);
    accel2(:,i) = myHighGrav.compAccelPar(reci(:,i)*1e3, 0);
err_mice(:,ii) = X_ecef(:,ii) - [recef_val(:,ii);vecef_val(:,ii)];
end
err_acc = accel-accel2;
%% Results: errors in in the order 2-3 m