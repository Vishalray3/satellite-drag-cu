%% Initialization

%% Furnishing spice kernels to use the spice functions
cspice_furnsh('de430.bsp')                       %% Planetary ephemerides
cspice_furnsh('naif0012.tls')                    %% time leap seconds
cspice_furnsh('earth_assoc_itrf93.tf')           %% Earth ITRF frame
cspice_furnsh('pck00010.tpc')
cspice_furnsh('gm_de431.tpc')                    %% GM values
cspice_furnsh('earth_000101_181111_180820.bpc')
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
epoch_utc     = 'Mar 18, 2017 10:29:00 UTC';            %% current epoch
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
doy = 178;                                     %% day of year (nrlmsise-00)
year = 2017;                                %% year (nrlmsise00)

eps  = 1000;                                   %% altitude diff. in m

%% Initial conditions
Hp = 300*1e3;            % perigee altitude
Ha = 500*1e3;
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
str_special = 'NO'; % CI: Circular inclined ; CE: circular equatorial; EE: elliptical equatorial; NO: None
[r_eci, v_eci] = coe2rv(a_sma,e,inc,raan,w_arg,true_ano,mu_e,'NO', u_arg);
n_circ = v_circ/r_circ;
T_orb = 2*pi/n_circ;
X_init = [r_eci;v_eci];

X_ref_st(:,1) = X_init + del_X(:,1);

%% Density model
switch flag_rho
    case 'Exp'
        g_acc = mu_e/r_p^2;          % gravity acceleration
        k_b = 1.38064852e-23; % m2 kg s-2 K-1
        amu = 1.66053904e-27; % kg
        atm_mass = [4,16,28,32,39,1,14]; % amu's
        
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
        load NRLMSISE_2007
        parameters.Ap_total = Ap_total;
        parameters.Ap_daily = Ap_daily;
        parameters.F10_total = F10_total;
    case 'JB08'
        load SOLFSMY_2019
        load DTCFILE_2019
        shape_model = wgs84Ellipsoid;                  %% shape model for geodetic to geocentric latitude
        flattening = shape_model.Flattening;           %% flattening factor
        ind_sol = find(SOLdata(1,:) == year,1);        %% index to point towards data for doy
        ind_mag = find(DTCdata(1,:) == year,1);                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        parameters.flattening = flattening;
        parameters.ind_sol = ind_sol;
        parameters.ind_mag = ind_mag;
        parameters.SOLdata = SOLdata;
        parameters.DTCdata = DTCdata;
end
%% SRP Parameters
if strcmp(flag_srp, 'Cball') == 1
    Cr = 1.9;
elseif strcmp(flag_srp, 'Three') == 1
    Cr = 1;
    A0 = -4.0;
    A1 = -1.6;
    A2 = -0.5;
end

%% Nominal values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch case_run
    case 'truth'
        if strcmp(flag_drag, 'Orb')
            Af_total = Af_total_mat(:,K_ind)';
            Af_true = Af_total;
            Bf_total = zeros(Order_b+1,Order_o+1);
            Bf_true = Bf_total;
        elseif strcmp(flag_drag, 'Bod')
            Af_total = Af_total_mat(:,:,K_ind);
            Af_true = Af_Kl(:,K_ind);           % Analytically averaged Cd
            Bf_total = Bf_total_mat(:,:,K_ind);
            Bf_true = Bf_Kl(:,K_ind);
        end
        Cf_total = zeros(Order_b+1,Order_o+1);
        Df_total = zeros(Order_b+1,Order_o+1);
        Cf_true = Cf_total; Df_true = Df_total;
    case {'estimation','consider'}
        load(truth_model,'X_true','R_aug','Af_true','Bf_true','Cf_true','Df_true','Cd_true')
        
        % Monte Carlo generated measurements
     for ii = 1:numel(time_pred)
        TI2E = cspice_sxform(  'J2000','ITRF93', time_pred(ii) );                      %% ECEF to ECI rot matrix for pos and vel
        GPS_state = TI2E*X_true(1:6,ii);
        GPS_pos(:,ii) = GPS_state(1:3) + Pos_noise(:,ii);
        GPS_vel(:,ii) = GPS_state(4:6) + Vel_noise(:,ii);
        time(ii) = time_pred(ii);

    end
    % Which measurement to use
    meas = 1:6;
    Measurements = [GPS_pos; GPS_vel];
    y_meas = Measurements(meas,:); 
        Af_total = Af_true + Af_noise ; % 
%         Af_total(:,1) =  Af_true(:,1)  + Af_noise(:,1); %
        Bf_total = Bf_true + Bf_noise ;
        Cf_total = Cf_true + Cf_noise ;
        Df_total = Df_true + Df_noise ;
end

X_f = [];
Xf_std = [];
Xf_true = [];
Cd = Af_total(1,1) ;%
Af0_true = Af_true(1,1);
if ismember('Cd',estimated_coeff)
    Cd_std = Af_std(1);
    % non-zero indices within given orders
    Af_ind = ~~Af_total(1:Order_b+1,1:Order_o+1);
    Bf_ind = ~~Bf_total(1:Order_b+1,1:Order_o+1);
    Cf_ind = ~~Cf_total(1:Order_b+1,1:Order_o+1);
    Df_ind = ~~Df_total(1:Order_b+1,1:Order_o+1);
    % nominal entries corresponing to the indices
    X_f = [Af_total(Af_ind); Bf_total(Bf_ind); Cf_total(Cf_ind); Df_total(Df_ind)];
    % nominal standard deviations
    Xf_std = [Af_std(Af_ind);  Bf_std(Bf_ind);  Cf_std(Cf_ind);  Df_std(Df_ind)];
    % truth values
    Xf_true = [Af_true(Af_ind); Bf_true(Bf_ind); Cf_true(Cf_ind); Df_true(Df_ind)];
    % Names of the coefficients
    mat_b = repmat([0:Order_b]',1,Order_o+1); mat_o = repmat([0:Order_o],Order_b+1,1);
    Af_name = strcat('A', sprintfc('%d',mat_b),sprintfc('%d',mat_o));
    Bf_name = strcat('B', sprintfc('%d',mat_b),sprintfc('%d',mat_o));
    Cf_name = strcat('C', sprintfc('%d',mat_b),sprintfc('%d',mat_o));
    Df_name = strcat('D', sprintfc('%d',mat_b),sprintfc('%d',mat_o));
    Xf_name_all = [Af_name(Af_ind); Bf_name(Bf_ind); Cf_name(Cf_ind); Df_name(Df_ind)];
    % vectors from those entries, : operator ensures column vector
    X_f = X_f(2:end)';
    Xf_std = Xf_std(2:end)';
    Xf_true = Xf_true(2:end)';
    Xf_name = Xf_name_all(2:end)';    
end

ind_est = find(ismember(Xf_name_all,str_est));
vec_est = [1:6,ind_est'+6];
if strcmp(str_case,'ignore')
X_f(setdiff(1:numel(X_f),vec_est(8:end)-7)) = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember('Cr',estimated_coeff) && ismember('Cd',estimated_coeff)
    if strcmp(flag_srp, 'Cball')
        N_nom = 8;
        X_nom = [X_ref_st;Cr;Cd;X_f'];
        P_prior = diag([del_r',del_v',Cr_std,Cd_std,Xf_std].^2);
    elseif strcmp(flag_srp, 'Three')
        X_nom = [X_ref_st;A0;A1;A2;Cd;X_f'];
        N_nom = 10;
        P_prior = diag([del_r',del_v',Cr_std,Cd_std,Xf_std].^2);
    end
elseif ismember('Cd',estimated_coeff)
    N_nom = 7;
    X_nom = [X_ref_st;Cd;X_f'];
    P_prior = diag([del_r',del_v',Cd_std,Xf_std].^2);
elseif isempty(estimated_coeff)
    N_nom = 6;
    X_nom = [X_ref_st];
    P_prior = diag([del_r',del_v'].^2);
end

N_st = N_nom + numel(X_f);
stm_init = eye(N_st);
x_prior = zeros(N_st,1);                 % initial deviation estimate
%% Force models
myTwoBody        = TwoBody;
myTwoBody.mu_e     = mu_e;

myDrag_atm         = Drag_atm;
myDrag_atm.Mass    = mass;
myDrag_atm.Area    = area;
myDrag_atm.flag    = flag_drag;
myDrag_atm.Ophi    = Order_o;
myDrag_atm.Cdrho   = 0;


Force(1).Model = myTwoBody;
Force(2).Model = myDrag_atm;

%% ODE options
ode4_delT = 10;
Tsamp = del_T/ode4_delT;
options = odeset('RelTol',1e-13,'AbsTol',1e-14);
step = time_pred_utc;
ode_step = step(1):ode4_delT:step(end);
%%
tym = step;
parameters.N_st = N_st;
parameters.N_nom = N_nom;
parameters.Order_o = Order_o;
parameters.Order_b = Order_b;
parameters.Af_ind = Af_ind;
parameters.Bf_ind = Bf_ind;
parameters.Cf_ind = Cf_ind;
parameters.Df_ind = Df_ind;
parameters.time_prop = time_prop;
parameters.tym = tym;
parameters.omega_e = omega_e;
parameters.epoch = epoch;
parameters.doy = doy;
parameters.year = year;
parameters.eps = eps;
parameters.Force = Force;
parameters.flag_srp = flag_srp;
parameters.flag_drag = flag_drag;
parameters.sun_pos = sun_pos;
parameters.moon_pos = moon_pos;
parameters.earth_vel = earth_vel;
parameters.bo_mod = bo_mod;
parameters.delT = del_T;
parameters.estimated_coeff = estimated_coeff;
parameters.flag_rho = flag_rho;
%% Truth or estimation run
if strcmp(case_run,'truth')
    xt = X_nom;
    xt(1:6,1) = X_init(:,1);
    X_true_aug(:,1) = [xt; stm_init(:)];
    X1 = ode4(@propagator_num,ode_step,X_true_aug,parameters,Tsamp);
    X_true = X1(:,1:6)';
    X_true_full = X1';
    
    for ii = 1:numel(time_pred)
        TI2E = cspice_sxform(  'J2000','ITRF93', time_pred(ii) );                      %% ECEF to ECI rot matrix for pos and vel
        
        Pos_noise(:,ii) = [normrnd(0,sigma_pos);normrnd(0,sigma_pos);normrnd(0,sigma_pos)];
        Vel_noise(:,ii) = [normrnd(0,sigma_vel);normrnd(0,sigma_vel);normrnd(0,sigma_vel)];
        GPS_state = TI2E*X_true(1:6,ii);
        GPS_pos(:,ii) = GPS_state(1:3) + Pos_noise(:,ii);
        GPS_vel(:,ii) = GPS_state(4:6) + Vel_noise(:,ii);
        time(ii) = time_pred(ii);
        
        [~,Cd_true(ii),rho_true(ii),body_angle(ii)] = propagator_num(time_pred_utc(ii),X_true_full(:,ii),parameters);
    end
    % Which measurement to use
    meas = 1:6;
    Measurements = [GPS_pos; GPS_vel];
    y_meas = Measurements(meas,:);
    R_aug = diag([sigma_pos,sigma_pos,sigma_pos,sigma_vel,sigma_vel,sigma_vel].^2);
    %     save('Truth_sphere_dria_Klast')
elseif strcmp(case_run,'estimation')
    run batch_obs
    Af_err_init = [Af0_true,Xf_true] - [Cd,X_f];
    Af_err_est = [Af0_true,Xf_true] - X_nom(7:end)';
    K_true = Kl_mat(K_ind);
    %     save(strcat(flag_drag, num2str(Ophi),'Area',num2str(area)))
     for ii = 1:numel(time_pred)
        [~,Cd_est(ii),~,~] = propagator_num_new(time_pred_utc(ii),X_aug(:,ii),parameters);
    end
    Cd_err = Cd_true - Cd_est;
    Cd_err_rms = rms(Cd_err);   
elseif strcmp(case_run,'consider')
    Af_err_init = [Af0_true,Xf_true] - [Cd,X_f];
    run batch_montecarlo
    Af_err_est = [Af0_true,Xf_true] - X_nom(7:end)';
    Xf_true_all = [Af0_true,Xf_true]; 
%     if strcmp(flag_drag,'Orb')
%     Af_err_est = Af_true(vec_est(7:end)-6) - X_nom(vec_est(7:end))';
%     elseif strcmp(flag_drag,'Bod_orb') || strcmp(flag_drag,'Bod')
%     Af_err_est = Xf_true(vec_est(7:end)-6) - X_nom(vec_est(7:end));
%     end
    K_true = Kl_mat(K_ind);
    Xstate_err = X_nom(1:6) - X_init;
    pos_err = norm(Xstate_err(1:3));
    vel_err = norm(Xstate_err(4:6));
    for ii = 1:numel(time_pred)
        [~,Cd_est(ii),~,~] = propagator_num(time_pred_utc(ii),X_aug(:,ii),parameters);
    end
    Cd_err = Cd_true - Cd_est;
    Cd_err_rms = rms(Cd_err);
    %     save(strcat(flag_drag, num2str(Ophi),'Area',num2str(area)))
end