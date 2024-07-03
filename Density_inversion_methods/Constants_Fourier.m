%% Constants for calculating the Fourier coefficients
%% Shape model
shape_model = 'sphere';           %% 'plate' or 'sphere'              
%% Furnishing spice kernels to use the spice functions
cspice_furnsh('de430.bsp')                       %% Planetary ephemerides
cspice_furnsh('naif0012.tls')                    %% time leap seconds
cspice_furnsh('earth_assoc_itrf93.tf')           %% Earth ITRF frame
cspice_furnsh('pck00010.tpc')
cspice_furnsh('gm_de431.tpc')                    %% GM values
cspice_furnsh('earth_000101_181111_180820.bpc')
%% Planetary constants
mu_e    = cspice_bodvrd('EARTH', 'GM',3)*1e9;           %% mu of Earth in m3/s2
Rad_e   = cspice_bodvrd('EARTH', 'RADII',3);            %% radius of Earth from the geophysical kernel (need for shadow function of SRP)
Re      = Rad_e(1)*1e3;
omega_e = 7.2921150e-5;                                 %% angular velocity of Earth in rad/s
%% Constants (Schamberg)
M_s = 27;                            % aluminum as surface mass
Ms = mean(M_s);
phi0 = 0;                             % Half angular width of beam of remitted molecules (0,10,25,45,60,90)
phi_s = goodman(phi0);                % gives the value of phi from lookup table defined by Goodman
%% Constants
R = 8314.46 ;                         % Universal gas constant (SI)
Tw = 300;                             % Temperature of the satellite wall
Kl = 1.44e6;                          % Mehta and walker .  % . Sesam 5e6/133.322;
Alpha = 1;
area_vec(:,1) = [1;0;0];

Area_plates = 1;      % solar panel areas = 13.5
Ar = 1;

phi = 0;
theta = 0;
flag_axis = 'z';
%% Fourier coefficients
N = 30;                                % highest order of coefficients
order_vec = [0:N]';
order_mat = repmat(order_vec,1,numel(theta));
theta_mat = order_mat.*theta;

%% Orbital information
case_ecc = 'low';
switch case_ecc
    case 'low'
        Hp = 300*1e3;            % perigee altitude
        Ha = 500*1e3;
    case 'high'
        Hp = 300*1e3;            % perigee altitude
        Ha = 36000*1e3;
end
r_p = Re+Hp;            % perigee radius
r_a = Re+Ha;
e = (r_a-r_p)/(r_a+r_p);   % eccentricity
a_sma = r_p/(1-e);        % semi major axis from perigee altitude
inc = 65;                  %%%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 97.8
raan = 60;
w_arg = 40;         % random
n_mean = sqrt(mu_e/a_sma^3);
true_ano = 0;      % random
u_arg = 0;          % argument of latitude
r_circ = a_sma;
v_circ = sqrt(mu_e/r_circ);
str_special = 'NO'; % CI: Circular inclined ; CE: circular equatorial; EE: elliptical equatorial; NO: None
[r_eci, v_eci] = coe2rv(a_sma,e,inc,raan,w_arg,true_ano,mu_e,'NO', u_arg);
n_circ = v_circ/r_circ;
T_orb = 2*pi/n_circ;
Fwind = (1-r_p*omega_e/norm(v_eci)*cosd(inc))^2;
%% Density calculations
global oxa amu k_b
g_acc = mu_e/r_p^2;          % gravity acceleration
k_b = 1.38064852e-23; % m2 kg s-2 K-1
amu = 1.66053904e-27; % kg
atm_mass = [4,16,28,32,39,1,14,16]; % amu's or molar mass:g/mol
he  = atm_mass(1);
oxa = atm_mass(2);
nitm= atm_mass(3);
oxm = atm_mass(4);
ar  = atm_mass(5);
h   = atm_mass(6);
nita= atm_mass(7);

[T_atm,n_den] = atmosnrlmsise00(Hp, 0, 0, 2019, 1, 0, 150, 150, 4*ones(1,7),'Oxygen');
% [T_atm, n_den] = atmosnrlmsise00(Hp, 0, 0, 2019, 1, 0,75, 75, 4*ones(1,7),'Oxygen');
% atmosnrlmsise00(altitude, latitude, longitude, year, dayOfYear, UTseconds, localApparentSolarTime);
rho0     = n_den(:,6);
rho0_he     = n_den(:,1)*he*amu;
rho0_oxa     = n_den(:,2)*oxa*amu;
rho0_nitm     = n_den(:,3)*nitm*amu;
rho0_oxm     = n_den(:,4)*oxm*amu;
rho0_ar     = n_den(:,5)*ar*amu;
rho0_h     = n_den(:,7)*h*amu;
rho0_nita     = n_den(:,8)*nita*amu;
rho0_oxa_ano = n_den(:,9)*oxa*amu;
rho0_all = [rho0_he,rho0_oxa,rho0_nitm,rho0_oxm,rho0_ar,rho0_h,rho0_nita,rho0_oxa_ano];
Talt = T_atm(2);

mass_sum = n_den(:,1)*he + n_den(:,2)*oxa + n_den(:,3)*nitm + n_den(:,4)*oxm + n_den(:,5)*ar + n_den(:,7)*h + n_den(:,8)*nita + n_den(:,9)*oxa;
num = sum(n_den,2);
M_mean = mass_sum./num*amu;

H_scale = k_b*Talt/(M_mean*g_acc);
H_he = k_b*Talt/(he*amu*g_acc);
H_oxa = k_b*Talt/(oxa*amu*g_acc);
H_nitm = k_b*Talt/(nitm*amu*g_acc);
H_oxm = k_b*Talt/(oxm*amu*g_acc);
H_ar = k_b*Talt/(ar*amu*g_acc);
H_h = k_b*Talt/(h*amu*g_acc);
H_nita = k_b*Talt/(nita*amu*g_acc);
H_scale_all = [H_he,H_oxa,H_nitm,H_oxm,H_ar,H_h,H_nita,H_oxa];