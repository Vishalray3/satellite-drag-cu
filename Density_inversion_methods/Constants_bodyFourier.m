%% Constants for calculating the Fourier coefficients
%% Furnishing spice kernels to use the spice functions
cspice_furnsh('de430.bsp')                       %% Planetary ephemerides
cspice_furnsh('naif0012.tls')                    %% time leap seconds
cspice_furnsh('earth_assoc_itrf93.tf')           %% Earth ITRF frame
cspice_furnsh('pck00010.tpc')
cspice_furnsh('gm_de431.tpc')                    %% GM values
cspice_furnsh('earth_000101_181111_180820.bpc')

% load('Cd_spire_fourier','Af_total_mat','Af_std','Bf_total_mat','Bf_std','Cf_total_mat','Cf_std','Df_total_mat','Df_std')
%% flags
flag_srp = 'Cball';                % Cball - cannonball, Three = three constants
if strcmp(case_run, 'truth')
    estimated_coeff = [{'CdTrue'}]; %[{'rho_DMC'},{'Cd'}];        % 'Cr','Cd','rho_DMC', 'cd_DMC', 'CdTrue', 'FourierDMC'
    flag_rho = 'MSIS00';                  % 'MSIS00', 'JB08', 'Exp'
    Kl = 1.44e6;    %1.44e6                      % Mehta and walker .  % . Sesam 5e6/133.322;
else
    estimated_coeff = [{'rho_DMC'},{'Cd'}];
    flag_rho = 'JB08';                  % 'MSIS00', 'JB08', 'Exp'
    Kl = 1e8;    %1.44e6                      % Mehta and walker .  % . Sesam 5e6/133.322;
end
flag_drag = 'Bod';                 % body model: Bod, orbit model: Orb; Body orbit model: Bod_orb
bo_mod = 'BOS';                   % BODF : don't forget to change both otheta and ophi
Order_b = 30;
Order_o = 0;
%% Planetary constants
mu_e    = cspice_bodvrd('EARTH', 'GM',3)*1e9;           %% mu of Earth in m3/s2
mu_s    = cspice_bodvrd('SUN', 'GM',3)*1e9;             %% mu of Sun in m3/s2
mu_m    = cspice_bodvrd('MOON', 'GM',3)*1e9;            %% mu of Moon in m3/s2
Rad_e   = cspice_bodvrd('EARTH', 'RADII',3);            %% radius of Earth from the geophysical kernel (need for shadow function of SRP)
Rad_s   = cspice_bodvrd('SUN', 'RADII',3);              %% radius of Earth from the geophysical kernel (need for shadow function of SRP)
Re      = Rad_e(1)*1e3;
Rs      = Rad_s(1)*1e3;
Psun    = 4.56e-6;                                      %% Radiation pressure at 1 au (Doornbos PhD) N/m2
AU      = 149597870660;                                 %% AU value m
omega_e = 7.2921150e-5;                                 %% angular velocity of Earth in rad/s
GPS_sec_offset = 19;                                    %% TAI to GPST offset (TAI-GPST)
deg_grav = 10 ;                                         %% degree of spherical harmonics
ord_grav = 10;                                          %% order of gravity harmonics
%% Time conversions  (ET epoch - Jan 1, 2000 12:00:00 TDB)
del_T         = 10;                                     %% sample time of measurements
tol_gap       = del_T;
epoch_utc     = 'Mar 23, 2007 00:00:00 UTC'; %'Oct 28, 2003 00:00:00 UTC'; %    %  %         %% current epoch
et_epoch      = cspice_str2et(epoch_utc);               %% convert current epoch in utc to et
tai_epoch     = cspice_unitim(et_epoch, 'TDB', 'TAI');  %% convert that to tai
time_prop_utc = 0:del_T:86400;  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 18
time_pred_utc = time_prop_utc;
tai_seconds   = tai_epoch + time_prop_utc;
tai_real      = tai_epoch + time_pred_utc;


time_prop     = cspice_unitim(tai_seconds, 'TAI', 'TDB');
time_pred     = cspice_unitim(tai_real, 'TAI', 'TDB');

%% Sun and Moon positions
[sun_pos, ~]  = cspice_spkpos('Sun', time_pred, 'J2000', 'NONE', 'EARTH');          %% positions of sun and moon in km
sun_pos = sun_pos*1e3;
[moon_pos, ~] = cspice_spkpos('Moon', time_pred, 'J2000', 'NONE', 'EARTH');
moon_pos = moon_pos*1e3;
[earth_state, ~] = cspice_spkezr('Earth', time_pred, 'J2000', 'NONE', 'Sun');
earth_state = earth_state*1e3;
earth_vel = earth_state(4:6,:);
%% GPM Parameters data for day 180 (June 29 2017)
mass = 500;                          %% in kg before the maneuver of Aug 9
area = 10;                                         %% change it for order 0
sun_area = 18.6;                                      %% reference area GPM, might need to change later
sun_area_mass = sun_area/mass;                        %% reference area for solar radiation  pressure

epoch = 0 ;                                    %% no. of seconds since the beginning of the day in UTC for the first observation (nrlmsise-00)
doy =   82;  % 302; %                     %% day of year (nrlmsise-00)
year = 2007;   %  2003; %                          %% year (nrlmsise00)

eps  = 1000;                                   %% altitude diff. in m

%% Body properties for SRP
R = 8314.46 ;                         % Universal gas constant (SI)
Tw = 300;                             % Temperature of the satellite wall
Alpha = 1;
area_vec(:,1) = [1;0;0];
area_vec(:,2) = [-1;0;0];
area_vec(:,3) = [0;1;0];
area_vec(:,4) = [0;-1;0];
area_vec(:,5) = [0;0;1];
area_vec(:,6) = [0;0;-1];
area_vec(:,7) = [0;0;1];               % feathered solar panel back side
area_vec(:,8) = [0;0;-1];              % feathered solar panel frontside
area_vec(:,9) = [0;0;1];               % feathered solar panel back side
area_vec(:,10) = [0;0;-1];             % feathered solar panel frontside

Area_plates = [1.5*2+1, 1.5*2+0.75, 2*2.5 + 1.5, 2*2.5 + 1.0, 2.5*1.5, 2.5*1.5+0.1, 13.5, 13.5, 13.5, 13.5]/area;      % solar panel areas = 13.5
% Area_plates = [1];      % sphere area
Ar = 1;

phi = 0;
theta = 0;
flag_axis = 'y';

M_s = [30,4,45,70,30,50,40,145,100,30];             % asymmetry test case
% M_s = [100];             % asymmetry test case
Ms = mean(M_s);
atm_mass = [4,16,28,32,39,1,14,16]; % amu's or molar mass:g/mol
amu = 1.66053904e-27; % kg
k_b = 1.38064852e-23; % m2 kg s-2 K-1
shape_model = 'plate_dria';  % sphere
%% Measurement noise values, actual GPM values
sig_meas(1) = 1.5; %3;  ; found by using the evar function on the residuals of the real data for order 150
sig_meas(2) = 5e-3; % 0.01; %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% noise mismatch
sigma_pos = sig_meas(1);
sigma_vel = sig_meas(2);

%% Initial conditions
Hp = 300*1e3;            % perigee altitude
Ha = 300*1e3;
% a_sma= Re+H;         % only for circular
r_p = Re+Hp;            % perigee radius
r_a = Re+Ha;
e = (r_a-r_p)/(r_a+r_p);   % eccentricity
a_sma = r_p/(1-e);        % semi major axis from perigee altitude
inc = 65;                  %%%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 97.8
raan = 0;
w_arg = 40;         % random
true_ano = 0;      % random
u_arg = 0;          % argument of latitude
r_circ = a_sma;
v_circ = sqrt(mu_e/r_circ);
str_special = 'CI'; % CI: Circular inclined ; CE: circular equatorial; EE: elliptical equatorial; NO: None
[r_eci, v_eci] = coe2rv(a_sma,e,inc,raan,w_arg,true_ano,mu_e,str_special, u_arg);
n_mean = v_circ/r_circ;
T_orb = 2*pi/n_mean;
Fwind = (1-r_p*omega_e/norm(v_eci)*cosd(inc))^2;
%% Fourier coefficients
order_vec = [0:Order_b]';
order_mat = repmat(order_vec,1,numel(theta));
theta_mat = order_mat.*theta;
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
        g_acc = mu_e/r_p^2;          % gravity acceleration
        [T_atm,n_den] = atmosnrlmsise00(Hp, 0, 0, 2019, 1, 0, 150, 150, 4*ones(1,7),'Oxygen');
        rho0     = n_den(:,6);
        Talt = T_atm(2);
        he  = atm_mass(1);                        %%%%% need to change to amu
        oxa = atm_mass(2);
        nitm= atm_mass(3);
        oxm = atm_mass(4);
        ar  = atm_mass(5);
        h   = atm_mass(6);
        nita= atm_mass(7);
        
        mass_sum = n_den(:,1)*he + n_den(:,2)*oxa + n_den(:,3)*nitm + n_den(:,4)*oxm + n_den(:,5)*ar + n_den(:,7)*h + n_den(:,8)*nita + n_den(:,9)*oxa;
        num = sum(n_den,2);
        M_mean = mass_sum./num*amu;
        
        H_scale = k_b*Talt/(M_mean*g_acc);
        parameters.rho0 = rho0;
        parameters.H_scale = H_scale;
        parameters.r0 = r_p;
    case 'MSIS00'
        load(strcat('NRLMSISE_',num2str(year)))
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