%% Load data

% load EIGEN-2
load GPS_data_2017_180
load Quat_eci2b_2017_180
load ephemeris_june27
load GravCoeff_EGM2008_zeroTide
% Cbar(3,1) = Cbar(3,1) -1.39e-8;
% load GPM_fourier_SP_canted
load('Cd_GPMreal_data180_fourier','Af_total_mat','Af_std','Bf_total_mat','Bf_std','Cf_total_mat','Cf_std','Df_total_mat','Df_std')
load('1_GSIM_MSIS_order2_iter1','Af_total_mat','Bf_total_mat')
% Af_total_mat(1,1) = 22.67;
% Af_total_mat(3,1) = -2.7;
% Af_total_mat(5,1) = 0.8;
% Af_total_mat(7,1) = -1.9;
% Af_std = Af_std/3;
% Bf_std = Bf_std/10;
%% Furnishing spice kernels to use the spice functions
cspice_furnsh('de430.bsp')                       %% Planetary ephemerides
cspice_furnsh('naif0012.tls')                    %% time leap seconds
cspice_furnsh('earth_assoc_itrf93.tf')           %% Earth ITRF frame
cspice_furnsh('pck00010.tpc')
cspice_furnsh('gm_de431.tpc')                    %% GM values
cspice_furnsh('earth_000101_210629_210407.bpc')

%% flags
vec_est = [1:36];
vec_meas = [1:6];                  % which measurements to consider
flag_srp = 'Cball';                % Cball - cannonball, Three = three constants, 'plate_model'
case_run = 'estimation';           % 'truth', 'estimation'
if strcmp(case_run, 'truth')
    estimated_coeff = [{'CdTrue'}]; %[{'rho_DMC'},{'Cd'}];        % 'Cr','Cd','rho_DMC', 'cd_DMC', 'CdTrue', 'FourierDMC'
    flag_rho = 'Exp';                  % 'MSIS00', 'JB08', 'Exp'
    Kl = 5e7;    %1.44e6                      % Mehta and walker .  % . Sesam 5e6/133.322;
else
    estimated_coeff = [{'Cr'},{'Cd'}]; %[{'rho_DMC'},{'Cd'}];
    flag_rho = 'MSIS00';                  % 'MSIS00', 'JB08', 'Exp'
    Kl = 5e7;    %1.44e6                      % Mehta and walker .  % . Sesam 5e6/133.322;
end
flag_drag = 'Bod';                 % body model: Bod, orbit model: Orb; Body orbit model: Bod_orb
bo_mod = 'BOS';                   % BODF : don't forget to change both otheta and ophi
Order_b = 30;
Order_o = 0;
Cd_nom = 0; %Af_total_mat(1,1);
flag_tides.SolidEarthTides = 1;
flag_tides.OceanTides = 1;
%% Planetary constants
c_light = 299792458.000000;                             %% Speed of light  [m/s]; DE430
mu_e    = cspice_bodvrd('EARTH', 'GM',3)*1e9;           %% mu of Earth in m3/s2 % 0.3986004415E+15; %
mu_s    = cspice_bodvrd('SUN', 'GM',3)*1e9;             %% mu of Sun in m3/s2
mu_m    = cspice_bodvrd('MOON', 'GM',3)*1e9;            %% mu of Moon in m3/s2
Rad_e   = cspice_bodvrd('EARTH', 'RADII',3);            %% radius of Earth from the geophysical kernel (need for shadow function of SRP)
Rad_s   = cspice_bodvrd('SUN', 'RADII',3);              %% radius of Earth from the geophysical kernel (need for shadow function of SRP)
Re      = Rad_e(1)*1e3; %0.6378136460E+07; %Rad_e(1)*1e3;
Rs      = Rad_s(1)*1e3;
Psun    = 4.56e-6;                                      %% Radiation pressure at 1 au (Doornbos PhD) N/m2
AU      = 149597870660;                                 %% AU value m
omega_e = 7.2921150e-5;                                 %% angular velocity of Earth in rad/s
GPS_sec_offset = 19;                                    %% TAI to GPST offset (TAI-GPST)
deg_grav = 150 ;                                         %% degree of spherical harmonics
ord_grav = 150;                                          %% order of gravity harmonics
%% Time conversions  (ET epoch - Jan 1, 2000 12:00:00 TDB)
end_index     = 17280;                                  %% last index till which data is being processed
time_string   = 'Jan 6, 1980 00:00:00 UTC';             %% GPS epoch
et_gps_epoch  = cspice_str2et(time_string);             %% convert gps epoch in utc to et
tai_gps_epoch = cspice_unitim(et_gps_epoch, 'TDB', 'TAI'); %% convert that to tai
sec_in_week   = 86400*7;                                %% seconds in a gps week
tai_time_gps  = GPS_week(1:end_index).*sec_in_week + GPST(1:end_index);  %% the epoch is in TAI and already is different from GPS by 19 s, no need to add offset
tai_time_gps  = tai_time_gps + tai_gps_epoch ;
et_time_gps   = cspice_unitim(tai_time_gps', 'TAI', 'TDB');
time_prop     = et_time_gps;
time_prop_utc = time_prop - time_prop(1);
%% Sun and Moon positions
[sun_pos, ~]  = cspice_spkpos('Sun', et_time_gps, 'J2000', 'NONE', 'EARTH');          %% positions of sun and moon in km
sun_pos = sun_pos*1e3;
[moon_pos, ~] = cspice_spkpos('Moon', et_time_gps, 'J2000', 'NONE', 'EARTH');
moon_pos = moon_pos*1e3;
[earth_state, ~] = cspice_spkezr('Earth', et_time_gps, 'J2000', 'NONE', 'Sun');
earth_state = earth_state*1e3;
earth_vel = earth_state(4:6,:);
%% GPM Parameters data for day 180 (June 29 2017)
mass = 3493.3981418762;                           %% in kg before the maneuver of Aug 9
area = 1;                                         %% change it for order 0
sun_area = 18.6;                                      %% reference area GPM, might need to change later
sun_area_mass = sun_area/mass;                        %% reference area for solar radiation  pressure



% X_ref_prev = [-3090449.948; 315765.3956; -6024350.586; -3292.78513; -6791.056671; 1333.092455];
% estimated state on Jun 28, 23:59:58.009 (from the definitive ephemeris file)

epoch_prev = 'Jun 28, 2017 23:59:58.009 UTC';
et_prev = cspice_str2et(epoch_prev);
delt = et_time_gps(1) - et_prev;
DV = datevec(epoch_prev); 
mjd_utc0 = mjuliandate(DV(1),DV(2),DV(3),DV(4),DV(5),DV(6)+delt);
t_req = 86400 + delt;                          %% see the ephemeris data to change this
x_req = interp1(sec_eph, x_eph, t_req);
y_req = interp1(sec_eph, y_eph, t_req);
z_req = interp1(sec_eph, z_eph, t_req);
vx_req = interp1(sec_eph, vx_eph, t_req);
vy_req = interp1(sec_eph, vy_eph, t_req);
vz_req = interp1(sec_eph, vz_eph, t_req);

reci = [x_req; y_req; z_req;];
veci = [vx_req; vy_req; vz_req];

epoch = 0 ;                                    %% no. of seconds since the beginning of the day in UTC for the first observation (nrlmsise-00)
doy = 180;                                     %% day of year (nrlmsise-00)
year = Year(1);                                %% year (nrlmsise00)
days_prev = 366;
eps  = 1000;                                   %% altitude diff. in m

%% read Earth orientation parameters
fid = fopen('eop19620101.txt','r');
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s 
%  ----------------------------------------------------------------------------------------------------
eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
fclose(fid);

ut1_utc = eopdata(7,:);
x_pole = eopdata(5,:);
y_pole = eopdata(6,:);
index_dut = find(eopdata(4,:)== floor(mjd_utc0));
%% Initial conditions
coe = rv2coe_E(reci(:,1),veci(:,1),mu_e);
a_sma = coe(1); e = coe(2); inc = coe(3); raan = coe(4); w_arg = coe(5); true_ano = coe(6); u_arg = coe(7);
if e > 1e-6
    str_special = 'CI';
else
    str_special = 'NO';
end
r_p = a_sma*(1-e);
Fwind = (1-r_p*omega_e/norm(veci(:,1))*cosd(inc))^2;
r_circ = a_sma;
v_circ = sqrt(mu_e/r_circ);
n_mean = v_circ/r_circ;
T_orb = 2*pi/n_mean;
Hp = a_sma*(1-e)-Re;
%% DMC model parameters
X_s = [0;0;0];
Xs_std = [1, 1e-2, 1e-1]; %[1e-12, 1e-13,1e-12];
zeta =  0.9;%0.2862;%0.2019; % 0.01; %
omega = 0.00145; %2*2*pi/T_orb; %6.5e-4; %0.0018; %2*pi/T_orb; % 2*2*pi/T_orb;
tau_inv = 0; %1/100;
sig_dmc = sqrt(5.5e-22);
q_dmc = 1;
c1 = 1e-5;%1e-17; % 1e-17; %sqrt(4*omega^3*zeta*sig_dmc^2/q_dmc);
c2 = 1.5e-5; %1e-18;%1e-16; %1e-16;
c3 = 1e-5;%1e-18;

% cd dmc parameters
X_cd = 0;
Xcd_std = 0.5;
tau_inv_cd = 0; %T_orb; %1/100;
c_cd = 0; % 0.5; %0.5e-3; %1e-3;
%% SRP Parameters
if strcmp(flag_srp, 'Cball') == 1
    Cr = 2.7;
    Cr_std = 1;
elseif strcmp(flag_srp, 'Three') == 1
    Cr = 1;
    A0 = -4.0;
    A1 = -1.6;
    A2 = -0.5;
else
    Cr = 1;
end
%% Body properties for drag and SRP
R = 8314.47215 ;                         % Universal gas constant (SI)
Tw = 300;                             % Temperature of the satellite wall
Alpha = 1;
area_vec(:,1) = [1;0;0];
area_vec(:,2) = [-1;0;0];
area_vec(:,3) = [0;1;0];
% area_vec(:,4) = [0;-1;0];
area_vec(:,4) = [0;0;1];
area_vec(:,5) = [0;0;-1];
area_vec(:,6) = [0;-sind(52);-cosd(52)];              % -y canted solar panel back side
area_vec(:,7) = [0;sind(52);cosd(52)];               % -y canted solar panel frontside
area_vec(:,8) = [0;0;-1];               % +y feathered solar panel back side
area_vec(:,9) = [0;0;1];             % +y feathered solar panel frontside

Area_plates = [4.2*2.5, 4.2*3, 3.5*2.5 + 2.5*0.5, 5*4.2, 5.*4.3, 13.5, 13.5, 13.5, 13.5];      % solar panel areas = 13.5
% Area_plates = [2*2.5, 2*2.5, 3.5*2.5 + 2.5*0.5, 5*4.2, 5.*4.3, 13.5, 13.5, 13.5, 13.5];      % solar panel areas = 13.5
% Area_plates = [0, 0, 0,0, 0, 13.5, 13.5, 13.5, 13.5];

Ar = 1;

spec_ref = [0.2, 0.05, 0.05, 0.2, 0.2, 0.2, 0.2, 0.05, 0.05];  % aluminum, Sio2 (solar panel) and polyurethane white paint A276 (solar panel backside)
diff_ref = [0.4, 0.3, 0.3, 0.4, 0.4, 0.4, 0.4, 0.3, 0.3];
N_plates = numel(spec_ref);
flag_axis = 'y';

M_s = [26.981538, 26.981538, 26.981538, 26.981538, 26.981538, 60.0843, 312.32, 60.0843, 312.32];    % satellite surface mass
Ms = mean(M_s);
atm_mass = [4.002602,15.9994,28.0134,31.9988,39.948,1.0079,14.0067, 15.9994]; % amu's or molar mass:g/mol; He, O, N2, O2, Ar, H, N, O
amu = 1.66053904e-27; % kg
k_b = 1.38064852e-23; % m2 kg s-2 K-1
shape_model = 'plate_dria';  % sphere

%% Density calculations
g_acc = mu_e/r_p^2;          % gravity acceleration
he  = atm_mass(1);
oxa = atm_mass(2);
nitm= atm_mass(3);
oxm = atm_mass(4);
ar  = atm_mass(5);
h   = atm_mass(6);
nita= atm_mass(7);

%% Density model
switch flag_rho
    case 'Exp'
        load nrlmsis_data_2016-17
        g_acc = mu_e/r_p^2;          % gravity acceleration
        [T_atm,n_den] = atmosnrlmsise00(Hp, 0, 0, year, 1, 0, F10_total(doy), F10_total(doy), 4*ones(1,7),'Oxygen');
        rho0     = n_den(:,6);
        Talt = T_atm(2);
        he  = atm_mass(1);                        %%%%% need to change to amu
        oxa = atm_mass(2);
        nitm= atm_mass(3);
        oxm = atm_mass(4);
        ar  = atm_mass(5);
        h   = atm_mass(6);
        nita= atm_mass(7);
        rho0_he     = n_den(:,1)*he*amu;
        rho0_oxa     = n_den(:,2)*oxa*amu;
        rho0_nitm     = n_den(:,3)*nitm*amu;
        rho0_oxm     = n_den(:,4)*oxm*amu;
        rho0_ar     = n_den(:,5)*ar*amu;
        rho0_h     = n_den(:,7)*h*amu;
        rho0_nita     = n_den(:,8)*nita*amu;
        rho0_oxa_ano = n_den(:,9)*oxa*amu;
        rho0_all = [rho0_he,rho0_oxa,rho0_nitm,rho0_oxm,rho0_ar,rho0_h,rho0_nita,rho0_oxa_ano];
        
        H_he = k_b*Talt/(he*amu*g_acc);
        H_oxa = k_b*Talt/(oxa*amu*g_acc);
        H_nitm = k_b*Talt/(nitm*amu*g_acc);
        H_oxm = k_b*Talt/(oxm*amu*g_acc);
        H_ar = k_b*Talt/(ar*amu*g_acc);
        H_h = k_b*Talt/(h*amu*g_acc);
        H_nita = k_b*Talt/(nita*amu*g_acc);
        H_scale_all = [H_he,H_oxa,H_nitm,H_oxm,H_ar,H_h,H_nita,H_oxa];
        mass_sum = n_den(:,1)*he + n_den(:,2)*oxa + n_den(:,3)*nitm + n_den(:,4)*oxm + n_den(:,5)*ar + n_den(:,7)*h + n_den(:,8)*nita + n_den(:,9)*oxa;
        num = sum(n_den,2);
        M_mean = mass_sum./num*amu;
        
        H_scale = k_b*Talt/(M_mean*g_acc);
        parameters.rho0 = rho0;
        parameters.H_scale = H_scale;
        parameters.r0 = r_p;
        parameters.rho0_all = rho0_all;
        parameters.H_scale_all = H_scale_all;
        parameters.Talt = Talt;
    case 'MSIS00'
        %         load(strcat('NRLMSISE_',num2str(year)))
        load nrlmsis_data_2016-17
        parameters.Ap_total = Ap_total;
        parameters.Ap_daily = Ap_daily;
        parameters.F10_total = F10_total;
        parameters.days_year_prev = 0;
    case 'JB08'
        load SOLFSMY_2019
        load DTCFILE_2019
        earth_model = wgs84Ellipsoid;                  %% shape model for geodetic to geocentric latitude
        flattening = earth_model.Flattening;           %% flattening factor
        ind_sol = find(SOLdata(1,:) == year,1);        %% index to point towards data for doy
        ind_mag = find(DTCdata(1,:) == year,1);                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        parameters.flattening = flattening;
        parameters.ind_sol = ind_sol;
        parameters.ind_mag = ind_mag;
        parameters.SOLdata = SOLdata;
        parameters.DTCdata = DTCdata;
end

%% Measurement noise values, actual GPM values
sig_meas(1) = 1.5; %3;  ; found by using the evar function on the residuals of the real data for order 150
sig_meas(2) = 5e-3; % 0.01; %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% noise mismatch
sigma_pos = sig_meas(1);
sigma_vel = sig_meas(2);
R_aug = diag([sigma_pos,sigma_pos,sigma_pos,sigma_vel,sigma_vel,sigma_vel].^2);
GPS_to_CM_offset_b = -[-0.5;1.0;-2.5];    %% In body fixed frame
GPS_pos = [GPS_ECEF_X(1:end_index)';GPS_ECEF_Y(1:end_index)';GPS_ECEF_Z(1:end_index)'];                  %% GPS positions in ITRF frame m
GPS_vel = [GPS_ECEF_VX(1:end_index)';GPS_ECEF_VY(1:end_index)';GPS_ECEF_VZ(1:end_index)'];               %% GPS velocities in ITRF frame m/s
GPS_meas_ecef = [GPS_pos; GPS_vel];
%% Using the quaternion data
q_ecib = [ q_ecib4(1:5:end), q_ecib1(1:5:end), q_ecib2(1:5:end), q_ecib3(1:5:end)];            %% first component is scalar
dcm_ecib = quat2dcm(q_ecib);                               %% dcm from quaternions

for ii = 1:end_index
    TIE = cspice_sxform( 'ITRF93', 'J2000', time_prop(ii) );                      %% ECEF to ECI rot matrix for pos and vel
    GPS_to_CM_offset_eci = dcm_ecib(:,:,ii)'*GPS_to_CM_offset_b;                  %% offset vector in eci frame
    GPS_meas_eci(:,ii) = TIE*GPS_meas_ecef(:,ii);
    GPS_pos_corr_eci(:,ii) = GPS_meas_eci(1:3,ii) + GPS_to_CM_offset_eci;
    RIE = cspice_pxform( 'ITRF93', 'J2000', time_prop(ii) );                      %% ECEF to ECI rot matrix for pos
    GPS_pos_corr(:,ii) = RIE'*GPS_pos_corr_eci(:,ii);
    omega_vec = RIE*[0;0;1]*omega_e;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     t = step(ii);
    %     UTsec = mod(epoch+t,86400);
    %     n_days = floor((epoch+t)/86400);
    %     day = doy + n_days;
    %
    %     X_ecef = GPS_pos_corr(:,ii);
    %     X_lla = ecef2lla(X_ecef');
    %     latitude = X_lla(1);
    %     longitude = X_lla(2);
    %     altitude = X_lla(3);
    %
    %     F10d = F10_total(day);                                   %% see nrlmsise_inputs_processing.m daily index
    %     F10  = mean(F10_total((day-40):(day+40)));                                   %% average of 81 days centered on current day
    %
    %
    %     UThour  = UTsec/3600;
    %     ap_daily = Ap_daily(day);                                  %% see nrlmsise_inputs_processing.m  ap daily index
    %     ap_ind = (day-1)*8+ceil(UThour/3);
    %     ap      = [ap_daily, Ap_total(ap_ind), Ap_total(ap_ind-1), Ap_total(ap_ind-2), Ap_total(ap_ind-3), ...
    %         mean(Ap_total((ap_ind-11):(ap_ind-4))), mean(Ap_total((ap_ind-19):(ap_ind-12)))];
    %     ap(~any(ap,1)) = ap_daily;
    %     wind = atmoshwm(latitude,longitude,altitude,'day',day,'seconds',UTsec,'apindex',Ap_total(ap_ind),'model','total');
    %     R_enu2ecef = [-sind(longitude), -cosd(longitude)*sind(latitude), cosd(longitude)*cosd(latitude);...
    %         cosd(longitude), -sind(longitude)*sind(latitude), sind(longitude)*cosd(latitude);...
    %         0, cosd(latitude), sind(latitude)];
    %     wind_ecef = R_enu2ecef*[wind(2); wind(1); 0];
    %     wind_eci = RIE*wind_ecef;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V_eci = GPS_meas_eci(4:6,ii) - cross(omega_vec, GPS_pos_corr_eci(:,ii));
    V_b = dcm_ecib(:,:,ii)*V_eci;
    SBF2ECI(:,:,ii) = dcm_ecib(:,:,ii)';
    if ismember('CdTrue',estimated_coeff)
        theta(1:numel(Area_plates),ii) = atan2d(V_b(3),V_b(1));
        phi(1:numel(Area_plates),ii) = atan2d(V_b(2), sqrt(V_b(1)^2+V_b(3)^2));
        theta(6:7,ii) = theta(6:7,ii) - SP_MY(5*ii-4);            % -y-axis canted panel
        theta(8:9,ii) = theta(8:9,ii) + SP_PY(5*ii-4);            % +y-axis feathered panel
    else
        theta(ii) = atan2d(V_b(3),V_b(1)) - SP_MY(5*ii-4);
        phi = 0;
    end
    angle_panels(ii,1,:) = zeros(numel(Area_plates),1);       % angle of panels w.r.t the body frame for panel based SRP
    angle_panels(ii,1,6:7) = SP_MY(5*ii-4);
    angle_panels(ii,1,8:9) = SP_PY(5*ii-4);
end
reci = GPS_pos_corr_eci;
veci = GPS_meas_eci(4:6,:);
Measurements = [GPS_pos_corr; GPS_vel];
y_meas = Measurements(vec_meas,:);
R_aug = R_aug(vec_meas,vec_meas);

%% Initialize state
del_r = [1;1;1];     % initial standard deviations
del_v = [0.1;0.1;0.1];
% X_init = [reci(:,1);veci(:,1)];
X_init = [GPS_pos_corr_eci(:,1);GPS_meas_eci(4:6,1)];

X_ref_st(:,1) = X_init;
% y_meas = [reci;veci];     % measurements
% y_meas = y_meas(vec_meas,:);
%% Initialize constants
parameters.e = e;
parameters.a_sma = a_sma;
parameters.inc = inc;
parameters.raan = raan;
parameters.w_arg = w_arg;
parameters.true_ano = true_ano;
parameters.u_arg = u_arg;
parameters.str_special = str_special;
parameters.Fwind = Fwind;
parameters.Re = Re;
parameters.mu_e = mu_e;
parameters.R = R;
parameters.Kl = Kl;
parameters.atm_mass = atm_mass;
parameters.M_s = M_s;
parameters.Ar = Ar;
parameters.Area_plates = Area_plates;
parameters.Tw = Tw;
parameters.Alpha = Alpha;
parameters.Order_b = Order_b;
parameters.Order_o = Order_o;
parameters.flag_axis = flag_axis;
parameters.area_vec = area_vec;
parameters.amu = amu;
parameters.sun_pos = sun_pos;
parameters.epoch = epoch;
parameters.doy = doy;
parameters.year = year;
parameters.eps = eps;
parameters.time_et = time_prop(1);
parameters.n_mean = n_mean;
parameters.flag_rho = flag_rho;
parameters.k_b = k_b;
parameters.shape_model = shape_model;
parameters.omega_e = omega_e;
parameters.rot_SBF2ECI = SBF2ECI;
% parameters.rot_ECI2ECEF = rot_ECI2ECEF;
parameters.theta = theta;
% parameters.eqeterms = eqeterms;
% parameters.eop = eop;
parameters.angle_panels = angle_panels;
parameters.phi = phi;
parameters.Re = Re;
parameters.Mjd_UTC0 = mjd_utc0;
parameters.dut = ut1_utc;
parameters.index_dut = index_dut;
parameters.flag_tides = flag_tides;
parameters.mu_m = mu_m;
parameters.mu_s = mu_s;
parameters.x_pole = x_pole;
parameters.y_pole = y_pole;