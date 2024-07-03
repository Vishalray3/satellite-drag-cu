% clc
% clearvars
restoredefaultpath
flag_disp = 0;
% dataset_sat = 'simulation';
%% Add path to folders
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/DataTools/mice/mice/src/mice/')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/DataTools/mice/mice/lib/' )
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/JB08')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods/data/ancillary_data')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods/data/HASDM_data')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/PosVelTransformations')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods/data/erp_data')
% 
% addpath('/home/vira0155/DataTools/mice/mice/src/mice/')
% addpath('/home/vira0155/DataTools/mice/mice/lib/' )
% addpath('/home/vira0155/Main_working_folder/JB08')
% addpath('/home/vira0155/Main_working_folder/Density_inversion_methods')
% addpath('/home/vira0155/Main_working_folder/Density_inversion_methods/data/ancillary_data')
% addpath('/home/vira0155/Main_working_folder/Density_inversion_methods/data/HASDM_data')
% addpath('/home/vira0155/Main_working_folder/PosVelTransformations')
% addpath('/home/vira0155/Main_working_folder/Density_inversion_methods/data/erp_data')

%% Flags
if strcmp(dataset_sat, 'spire')
    Constants_spire
    propagator_dmc_real = @propagator_dmc_real_spire;
    ekf_script = 'ekf_dmc_smoothing_real';
elseif strcmp(dataset_sat, 'starlink')
    Constants_starlink
    propagator_dmc_real = @propagator_dmc_real_starlink;
    ekf_script = 'ekf_dmc_smoothing_real';
elseif strcmp(dataset_sat, 'simulation')
    Constants_simSpire
    propagator_dmc_real = @propagator_dmc;
    ekf_script = 'ekf_dmc_smoothing';
end
%% Nominal values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch case_run
    case 'Truth'
        if strcmp(flag_drag, 'Orb')
            Af_total = Af_total_mat(:,K_ind)';
            Af_true = Af_total;
            Bf_total = zeros(Order_b+1,Order_o+1);
            Bf_true = Bf_total;
            Cf_total = zeros(Order_b+1,Order_o+1);
            Cf_true = Cf_total;
            Df_total = zeros(Order_b+1,Order_o+1);
            Df_true = Df_total;
            Cf_std = zeros(Order_b+1,Order_o+1);
            Df_std = zeros(Order_b+1,Order_o+1);
        elseif strcmp(flag_drag, 'Bod') || strcmp(flag_drag, 'Bod_orb')
            Af_total = Af_total_mat(:,:,K_ind);
            Af_true  = Af_total; %Af_Kl(:,K_ind);           % Analytically averaged Cd
            Bf_total = Bf_total_mat(:,:,K_ind);
            Bf_true  = Bf_total; %Bf_Kl(:,K_ind);
            Cf_total = Cf_total_mat(:,:,K_ind);
            Df_total = Df_total_mat(:,:,K_ind);
            Cf_true = Cf_total; Df_true = Df_total;
        end
        
    case {'Estimation','consider'}
        %         load(truth_model,'X_true','R_aug','Cd_true','rho_true','Af_tseries','Bf_tseries')
        Af_total = Af_total_mat(:,:,K_ind);
        Bf_total = Bf_total_mat(:,:,K_ind);
        Cf_total = Cf_total_mat(:,:,K_ind);
        Df_total = Df_total_mat(:,:,K_ind);
        if gsim_iter == 2
            load(truth_model, 'Af_pred','Bf_pred')
            Af_total = Af_pred;
            Bf_total = Bf_pred;
        end
end

X_f = [];
Xf_std = [];
Xf_true = [];
Cd = 0; %Cd_nom; %Af_total(1,1);
Af_ind = []; Bf_ind = []; Cf_ind = []; Df_ind = [];

if estimated_coeff.Cd || estimated_coeff.CdDiscrete_est
    Cd_std = Xcd_std; %Af_std(1);
    % non-zero indices within given orders
    Af_ind = ~~Af_total(1:Order_b+1,1:Order_o+1);
    Bf_ind = ~~Bf_total(1:Order_b+1,1:Order_o+1);
    Cf_ind = ~~Cf_total(1:Order_b+1,1:Order_o+1);
    Df_ind = ~~Df_total(1:Order_b+1,1:Order_o+1);
    % nominal entries corresponing to the indices
    Xf1 = Af_total(Af_ind); Xf2 = Bf_total(Bf_ind); Xf3 = Cf_total(Cf_ind); Xf4 = Df_total(Df_ind);
    X_f = [Xf1(:); Xf2(:); Xf3(:); Xf4(:)];
    % nominal standard deviations
    Xf1 = Af_std(Af_ind); Xf2 = Bf_std(Bf_ind); Xf3 = Cf_std(Cf_ind); Xf4 = Df_std(Df_ind);
    Xf_std= [Xf1(:); Xf2(:); Xf3(:); Xf4(:)];
    % truth values
    %     Xf1 = Af_true(Af_ind); Xf2 = Bf_true(Bf_ind); Xf3 = Cf_true(Cf_ind); Xf4 = Df_true(Df_ind);
    %     Xf_true= [Xf1(:); Xf2(:); Xf3(:); Xf4(:)];
    % Names of the coefficients
    mat_b = repmat([0:Order_b]',1,Order_o+1); mat_o = repmat([0:Order_o],Order_b+1,1);
    Af_name = strcat('A', sprintfc('%d',mat_b),sprintfc('%d',mat_o));
    Bf_name = strcat('B', sprintfc('%d',mat_b),sprintfc('%d',mat_o));
    Cf_name = strcat('C', sprintfc('%d',mat_b),sprintfc('%d',mat_o));
    Df_name = strcat('D', sprintfc('%d',mat_b),sprintfc('%d',mat_o));
    Xf1_name = Af_name(Af_ind); Xf2_name = Bf_name(Bf_ind); Xf3_name = Cf_name(Cf_ind); Xf4_name = Df_name(Df_ind);
    Xf_name_all = [Xf1_name(:); Xf2_name(:); Xf3_name(:); Xf4_name(:)];
    % vectors from those entries, : operator ensures column vector
    
    X_f = X_f(2:end)';
    Xf_std = Xf_std(2:end)';
    %     Xf_true = Xf_true(2:end)';
    Xf_name = Xf_name_all(2:end)';
    
    Af_ind1 = ~~Af_std(1:Order_b+1,1:Order_o+1);
    Bf_ind1 = ~~Bf_std(1:Order_b+1,1:Order_o+1);
    Cf_ind1 = ~~Cf_std(1:Order_b+1,1:Order_o+1);
    Df_ind1 = ~~Df_std(1:Order_b+1,1:Order_o+1);
    Xf1_name = Af_name(Af_ind1); Xf2_name = Bf_name(Bf_ind1); Xf3_name = Cf_name(Cf_ind1); Xf4_name = Df_name(Df_ind1);
    Xf_name_std = [Xf1_name(:); Xf2_name(:); Xf3_name(:); Xf4_name(:)];
    
end


%%
N_nom = 6;
X_nom = X_ref_st;
P_prior = diag([del_r',del_v'].^2);
X_name_nom = [{'p1'},{'p2'},{'p3'},{'v1'},{'v2'},{'v3'}];
if estimated_coeff.Cr
    if strcmp(flag_srp, 'Cball') || strcmp(flag_srp, 'Panel')
        N_nom = N_nom+1;
        X_nom = [X_nom;Cr];
        P_prior = blkdiag(P_prior,Cr_std^2);
        X_name_nom = [X_name_nom,{'Cr'}];
    elseif strcmp(flag_srp, 'Three')
        X_nom = [X_ref_st;A0;A1;A2];
        N_nom = N_nom+3;
        P_prior = blkdiag(P_prior,diag([A0_std^2, A1_std^2, A2_std^2]));
        X_name_nom = [X_name_nom,{'A0'},{'A1'},{'A2'}];
    end
end

stm_q_init = [];
if estimated_coeff.rho_DMC
    N_nom = N_nom+3;
    X_nom = [X_nom; X_s];
    P_prior = blkdiag(P_prior, diag(Xs_std.^2));
    X_name_nom = [X_name_nom,{'rho1'},{'rho2'},{'rho3'}];
    Nrho = N_nom;
end

if estimated_coeff.Cd
    N_nom = N_nom+1;
    X_nom = [X_nom;Cd;X_f'];
    P_prior = blkdiag(P_prior, diag([Cd_std^2, Xf_std.^2]));
    X_name_nom = [X_name_nom,Xf_name_all'];
end
if estimated_coeff.CdDiscrete_est
    N_nom = N_nom+1;
    X_nom = [X_nom;Cd_nom;X_f'];
    P_prior = blkdiag(P_prior, diag([Cd_std^2, Xf_std.^2]));
end

if estimated_coeff.Cr_erp
    N_nom = N_nom+1;
    X_nom = [X_nom;Cr_erp];
    X_name_nom = [X_name_nom, {'Cr_erp'}];
    P_prior = blkdiag(P_prior, Cr_erp_std^2);
end

if parameters.empirical
    X_nom = [X_nom; X_emp];
    X_name_nom = [X_name_nom, {'accn'},{'acct'},{'accw'}];
    P_prior = blkdiag(P_prior, diag(Xemp_std.^2));
    N_emp = 3;
else
    N_emp = 0;
end

N_f = numel(X_f);
parameters.N_f = N_f;
N_st = N_nom + N_f + N_emp;
stm_init = eye(N_st);
if strcmp(case_run, 'Truth')
    stm_init = stm_init(vec_est,vec_est);
end
x_prior = zeros(N_st,1);                 % initial deviation estimate
%
% if exist('X_nom_new','var')
%     X_nom(9) = X_nom_new(9);
% end
if estimated_coeff.rho_DMC
    stm_q_init = zeros(N_st);
end
X_name_est = X_name_nom(vec_est);
%% ODE options
ode4_delT = del_T;
Tsamp = del_T/ode4_delT;
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
step = time_prop_utc;
ode_step = step(1):ode4_delT:step(end);
%%
tym = step;
parameters.N_st = N_st;
parameters.N_nom = N_nom;
parameters.Af_ind = Af_ind;
parameters.Bf_ind = Bf_ind;
parameters.Cf_ind = Cf_ind;
parameters.Df_ind = Df_ind;
parameters.time_prop = time_prop;
parameters.tym = tym;
parameters.epoch = epoch;
parameters.doy = doy;
parameters.year = year;
parameters.eps = eps;

parameters.sun_pos = sun_pos;
parameters.moon_pos = moon_pos;
parameters.earth_vel = earth_vel;
parameters.bo_mod = bo_mod;
parameters.delT = del_T;
parameters.estimated_coeff = estimated_coeff;
parameters.flag_rho = flag_rho;
parameters.vec_est = vec_est;
parameters.days_year_prev = days_prev;
parameters.zeta = zeta;
parameters.omega = omega;
parameters.tau_inv = tau_inv;
parameters.tau_inv_cd = tau_inv_cd;
parameters.c1 = c1;
parameters.c2 = c2;
parameters.c3 = c3;
parameters.c_cd = c_cd;
parameters.Cd_est = Cd_nom;
parameters.Cr_est = Cr;
parameters.T_orb = T_orb;
parameters.atm_mass = atm_mass;
parameters.amu = amu;
parameters.Cnm = Cbar(1:deg_grav+1, 1:deg_grav+1);
parameters.Snm = Sbar(1:deg_grav+1, 1:deg_grav+1);

% parameters.theta_max = theta_max;
if estimated_coeff.CdTrue
    parameters.R = R; parameters.Kl = Kl; parameters.M_s = M_s; parameters.Ar = Ar;
    parameters.Area_plates = Area_plates; parameters.Tw = Tw; parameters.Alpha = Alpha;  parameters.flag_axis = flag_axis;
    parameters.area_vec = area_vec; parameters.flag_rho = flag_rho; parameters.k_b = k_b; parameters.shape_model = shape_model;
end
%% Truth or estimation run
switch case_run
    case 'EDR'
        if strcmp(dataset_sat, 'spire') == 1
            run EDR_spire
        elseif strcmp(dataset_sat, 'simulation') == 1
            load(truth_model,'X_true_aug')
            reci = X_true_aug(1:3,:);
            veci = X_true_aug(4:6,:);
            run EDR_noise
        end
        
    case 'Truth'
        X_true_aug(:,1) = [X_init(:,1)];
        X1 = ode4(@propagator_truth,ode_step,X_true_aug,parameters,Tsamp);
        %     X_reg = X1(:,1:6)';
        %     X_true = interp1(ode_step, X_reg',time_pred_utc,'spline');
        %     X_true = X_true';
        %     X_true_full = X1';
        %         [t,X1] = ode45(@(t,x)propagator_truth(t, x, parameters), time_prop_utc, X_true_aug(:,1), options);
        X_true_aug = X1';
        %         deltaV_drag = X_true_aug(7, end);
        
        clear reci veci
        for ii = 1:numel(time_prop)
            %         TI2E = cspice_sxform(  'J2000','ITRF93', time_pred(ii) );                      %% ECEF to ECI rot matrix for pos and vel
            
            Pos_noise(:,ii) =[sigma_pos*randn;sigma_pos*randn;sigma_pos*randn];
            Vel_noise(:,ii) = [sigma_vel*randn;sigma_vel*randn;sigma_vel*randn];
            GPS_state = X_true_aug(1:6,ii);
            reci(:,ii) = GPS_state(1:3) + Pos_noise(:,ii);
            veci(:,ii) = GPS_state(4:6) + Vel_noise(:,ii);
            %             time(ii) = time_prop(ii);
            [~,Cd_est(ii),rho(ii), theta(ii), phi(ii),a_grav,a_sun,a_moon,a_drag(:,ii),a_srp(:,ii),a_earthrad(:,ii)] = ...
                propagator_truth(time_prop_utc(ii),X_true_aug(:,ii),parameters);
        end
        r_ind = round(T_orb/10);
        r_alt = vecnorm(reci(1:3,:),2,1)-Re;
        ii = 1;
        for kk=1:r_ind:numel(r_alt)-r_ind
            r_alt_avg(ii) = mean(r_alt(kk:kk+r_ind));
            r_alt_min(ii) = min(r_alt(kk:kk+r_ind));
            rho_avg(ii) = mean(rho(kk:kk+r_ind));
            ii = ii+1;
        end
        time_avg = time_prop_utc(r_ind:r_ind:end);
    case 'Parameters'
        X_full = zeros(2*N_st^2+2*N_st, numel(time_prop_utc_ekf));
        X_full(1:6,:) = [reci;veci];
        for ii = 1:numel(time_prop_utc_ekf)
            [~,Cd_est(ii),rho_est(ii),rho_nom(ii),a_srp(:,ii),a_earthrad(:,ii), a_drag(:,ii),Aref(ii)] = propagator_dmc_real(time_prop_utc_ekf(ii),X_full(:,ii),parameters);
        end
        Cd_ref = Cd_est./Aref;
        for ii = 1:numel(jdutc_sec)
            X_lla = ecef2lla(recef_ind(ii,:));
            latitude(ii) = X_lla(1);
            longitude(ii) = X_lla(2);
            height_calc(ii) = X_lla(3)/1e3;
        end
        clear year
        data_mat = [year(t_vec),month(t_vec),day(t_vec),hour(t_vec),minute(t_vec),second(t_vec),latitude', longitude', height_calc',Cd_ref', Aref',...
            recef_ind, vecef_ind, q_ecef2body, sada_angle];
        data_tab = array2table(data_mat);
        data_tab.Properties.VariableNames = {'Year','Month','Day','Hour(UTC)','Minute(UTC)','Second(UTC)','Latitude (deg)','Longitude (deg)',...
            'Altitude (km)','Cd','Aref (m2)', 'recef1 (m)','recef2 (m)','recef3 (m)','vecef1 (m/s)','vecef2 (m/s)','vecef3 (m/s)', 'qecef2body0 (scalar)'...
            , 'qecef2body1', 'qecef2body2', 'qecef2body3', 'solar panel angle (deg)'};
        writetable(data_tab, 'starlink3167_ephemeris_full.csv')
    case 'Estimation'
        if strcmp(case_est, 'EKF')
            
            if strcmp(dataset_sat, 'simulation')
%                 load(truth_model,'X_true_aug', 'rho')
                sigma_pos = sigma_pos_new;
                sigma_vel = sigma_vel_new;
                for ii = 1:numel(time_prop)
                    Pos_noise(:,ii) =[sigma_pos*randn;sigma_pos*randn;sigma_pos*randn];
                    Vel_noise(:,ii) = [sigma_vel*randn;sigma_vel*randn;sigma_vel*randn];
                    GPS_state = X_true_aug(1:6,ii);
                    reci(:,ii) = GPS_state(1:3) + Pos_noise(:,ii);
                    veci(:,ii) = GPS_state(4:6) + Vel_noise(:,ii);
                end
                y_meas = [reci; veci];
                y_meas = y_meas(vec_meas,:);
            end
            
            run(ekf_script)
            
            parameters_temp = parameters;
            parameters_temp.sun       = 0;
            parameters_temp.moon      = 0;
            parameters_temp.tides     = 0;
            parameters_temp.relativity= 0;
            flag_tides.SolidEarthTides = 0;
            flag_tides.OceanTides = 0;
            deg_grav = 0;
            parameters_temp.deg_grav = deg_grav;
            parameters_temp.ord_grav = deg_grav;
            parameters_temp.Cnm = Cbar(1:deg_grav+1, 1:deg_grav+1);
            parameters_temp.Snm = Sbar(1:deg_grav+1, 1:deg_grav+1);
            %             parameters.vec_est = [1:9];
            for ii = 1:numel(time_prop_utc_ekf)
                X_full(vec_est,ii) = Xs_est(:,ii);
                [~,Cd_est(ii),rho_est(ii),rho_nom(ii),a_total(ii), a_drag(:,ii), a_srp(:,ii),a_earth_rad(:,ii), rot_ECI2ECEF] = ...
                    propagator_dmc_real(time_prop_utc_ekf(ii),X_full(:,ii),parameters_temp);
                
                X_ecef = rot_ECI2ECEF*Xs_est(1:3,ii);
                X_lla = ecef2lla(X_ecef');
                latitude(ii) = X_lla(1);
                longitude(ii) = X_lla(2);
                alt_calc(ii) = X_lla(3);
%                 alt_geoid = egm96geoid(latitude(ii), longitude(ii));
%                 alt_calc(ii) = altitude_ellip - alt_geoid;     % geoidal height above EGM-96
                %                 sigma_rhoest(ii) = rho_nom(ii)*sqrt(sigma_x(Nrho-2,ii)^2 + sigma_x(Nrho,ii)^2 + 2*sigma_x(Nrho,ii)*sigma_x(Nrho-2,ii));
                %                 sigma_rhoests(ii) = rho_nom(ii)*sqrt(sigma_xs(Nrho-2,ii)^2 + sigma_xs(Nrho,ii)^2 + 2*sigma_xs(Nrho,ii)*sigma_xs(Nrho-2,ii));
            end
            
            if strcmp(dataset_sat, 'spire')
                timediff = diff(time_prop_utc_ekf);
                ind = find(timediff > 10);
                if ~isempty(ind)
                    ind_arc = [ind(1)+1:numel(rho_est)];
%                     ind_pos = find(rho_est>0);
%                     ind_vec = intersect(ind_arc, ind_pos);
                    ind_vec = ind_arc;
                    time_prop_trunc = time_prop_utc_ekf(ind_vec);
                    rho_est_trunc = rho_est(ind_vec);
                    rho_nom_trunc = rho_nom(ind_vec);
                    lat_trunc = latitude(ind_vec);
                    long_trunc = longitude(ind_vec);
                    alt_trunc = alt_calc(ind_vec);
                    jdutc_trunc = jdutc_sec_data(ind_vec);
                    X_state_trunc = Xs_est(1:6, ind_vec);
                    sun_pos_trunc = sun_pos_data(:,ind_vec);
                    moon_pos_trunc = moon_pos_data(:,ind_vec);
                    sod_trunc = sod_data(ind_vec);
                else
                    ind_vec = [];
                    time_prop_trunc = [];
                    rho_est_trunc = [];
                    rho_nom_trunc = [];
                    lat_trunc = [];
                    long_trunc = [];
                    alt_trunc = [];
                    jdutc_trunc = [];
                    X_state_trunc = [];
                    sun_pos_trunc = sun_pos_data(:,ind_vec);
                    moon_pos_trunc = moon_pos_data(:,ind_vec);
                    sod_trunc = [];
                end
                
%                 F_hasdm = hasdm_interpolant(n_days, doy, hasdm_mat);
%                 jd_ref = 86400*GREGORIANtoJD_vector(1950,1,0) + 3600*0 + 60*0 + 0;
%                 njd = (jdutc_trunc - jd_ref)/86400;
%                 long_trunc_pos = long_trunc;
%                 long_trunc_pos(long_trunc_pos<0) = long_trunc_pos(long_trunc_pos<0) + 360;
%                 rho_hasdm_trunc = exp(F_hasdm(alt_trunc/1e3, njd,long_trunc_pos, lat_trunc));
%                 save(strcat(dir_data,'/spireResults_', flag_rho, '_',sat_ID,'_', date, '_new'),'rho_est','Cd_est', 'rho_nom','x0_est_iter','y_meas',...
%                     'ys_res_postfit','time_prop_utc_ekf','time_prop_utc','X_nom', 'Xs_est', 'theta_calc','phi_calc', 'height_calc',...
%                     'a_sma', 'e', 'inc','raan', 'w_arg','yyyy','mon','day_data','hh','mm','ss','EOPInfo','NUT1980Info', 'LeapSecInfo','jdutc_sec_data', 'vec_est',...
%                     'estimated_coeff','latitude','longitude','alt_calc','ind_vec', 'time_prop_trunc','rho_est_trunc','rho_nom_trunc',...
%                     'lat_trunc','long_trunc','alt_trunc','jdutc_trunc',...
%                     'doy', 'sun_pos_trunc','moon_pos_trunc', 'epoch','X_state_trunc', 'rho_hasdm_trunc', 'a_earth_rad', 'a_drag', 'a_srp', 'sod_trunc')
            elseif strcmp(dataset_sat, 'starlink')
                
                ind_arc = [200:numel(rho_est)];
                ind_pos = find(rho_est>0);
                ind_vec = intersect(ind_arc, ind_pos);
                time_prop_trunc = time_prop_utc_ekf(ind_vec);
                rho_est_trunc = rho_est(ind_vec);
                rho_nom_trunc = rho_nom(ind_vec);
                lat_trunc = latitude(ind_vec);
                long_trunc = longitude(ind_vec);
                alt_trunc = alt_calc(ind_vec);
                jdutc_trunc = jdutc_sec(ind_vec);
                X_state_trunc = Xs_est(1:6, ind_vec);
                sun_pos_trunc = sun_pos(:,ind_vec);
                moon_pos_trunc = moon_pos(:,ind_vec);
                save(strcat(dir_data,'/starlinkResults_', flag_rho, '_',sat_ID,'_', date, '_new'),'rho_est','Cd_est', 'rho_nom','x0_est_iter','y_meas',...
                    'ys_res_postfit','time_prop_utc_ekf','time_prop_utc','X_nom', 'Xs_est', 'theta_calc','phi_calc',...
                    'a_sma', 'e', 'inc','raan', 'w_arg','yyyy','doy_vec','eop','jdutc_sec', 'vec_est',...
                    'estimated_coeff','latitude','longitude','alt_calc', 'time_prop','rho_est','rho_nom', 'Cd_s','Cd_ads', 'index_vec', ...
                    'ind_vec', 'time_prop_trunc','rho_est_trunc','rho_nom_trunc', 'lat_trunc','long_trunc','alt_trunc','jdutc_trunc',...
                    'doy', 'sun_pos_trunc','moon_pos_trunc', 'epoch','X_state_trunc')
                
            elseif strcmp(dataset_sat, 'simulation')
                rho_est_error = (rho_true - rho_est)./rho_true*100;
                rho_est_mean = mean(rho_est_error);
                rho_est_rms  = rms(rho_est_error);
                
                for ii = 1:numel(time_prop_utc_ekf)
                    
                    [~,~,~,~,~,a_grav,a_sun,a_moon,a_drag,a_srp,a_earthrad] = ...
                        propagator_truth(time_prop_utc(ii),X_true_aug(:,ii),parameters_truth);        
                    a_drag_true(ii) = norm(a_drag);              
                end
                [rho_est_mean_avg, rho_est_rms_avg] = den_avg(rho_true, rho_est, arclen);
                rms_min = min(rho_est_rms_avg);
                ind_min = find(rho_est_rms_avg < 1.1*rms_min);
                rho_rms_avg_min = rho_est_rms_avg(ind_min(1));
                rho_mean_avg_min = rho_est_mean_avg(ind_min(1));
                arc_min = arclen(ind_min(1));
                ind_10 = find(rho_est_rms_avg < 10);
                arc_10 = arclen(ind_10(1));
                ind_5 = find(rho_est_rms_avg < 5);
                arc_5 = arclen(ind_5(1));
                
                P_signal = rms(a_drag_true);
            end
        elseif strcmp(case_est, 'batch')
            run batch_considercovariance
        end
        
    case 'considercovariance'
        %     Af_err_init = [Af0_true,Xf_true] - [Cd,X_f];
        run batch_considercovariance
        %     Af_err_est = Af_true(vec_est(7:end)-6) - X_nom(vec_est(7:end))';
        K_true = Kl_mat(K_ind);
        Xstate_err = X_nom(1:6) - X_init;
        pos_err = norm(Xstate_err(1:3));
        vel_err = norm(Xstate_err(4:6));
        for ii = 1:numel(time_pred)
            [~,Cd_est(ii),~,~] = propagator_num(time_pred_utc(ii),X_aug(:,ii),parameters);
        end
        Cd_err = Cd_true - Cd_est;
        Cd_err_rms = rms(Cd_err);
end


%% Propagation
% % Change parameters
% scale_factor = mean(rho_est)/mean(rho_nom);
% Cd_nom = scale_factor*Cd_nom;
% parameters.Cd_est = Cd_nom;
% 
% 
% step = 0:10:86400*3;
% parameters.jd_init = jdutc_sec_data(end);
% tai_time_prop =  tai_time_gps_data(end) + step; 
% et_time_gps   = cspice_unitim(tai_time_prop, 'TAI', 'TDB');
% [sun_pos, ~]  = cspice_spkpos('Sun', et_time_gps, 'J2000', 'NONE', 'EARTH');          %% positions of sun and moon in km
% sun_pos = sun_pos*1e3;
% [moon_pos, ~] = cspice_spkpos('Moon', et_time_gps, 'J2000', 'NONE', 'EARTH');
% moon_pos = moon_pos*1e3;
% 
% [earth_state, ~] = cspice_spkezr('Earth', et_time_gps, 'J2000', 'NONE', 'Sun');
% earth_state = earth_state*1e3;
% earth_vel = earth_state(4:6,:);
% 
% parameters.sun_pos = sun_pos;
% parameters.moon_pos = moon_pos;
% parameters.earth_vel = earth_vel;
% parameters.time_et = et_time_gps;
% parameters.tym = step;
% 
% X_true_aug = X_full(:,1);
% X_true_aug(vec_est) = Xs_est(vec_est, end);
% 
% if parameters.estimated_coeff.rho_DMC 
%     X_true_aug(8:10) = 0;
% end
% parameters.vec_est = vec_est;
% 
% ode_step = step(1):ode4_delT:step(end);
% q_bo = zeros(4, numel(step));
% parameters.q_bo = q_bo;
% 
% 
% X1_prop = ode4(propagator_dmc_real,ode_step,X_true_aug(:,1),parameters,Tsamp);
% 
% 
% % Load truth
% date_curr = '2023_09_05';
% load(strcat('spire_satellite_data',data_pod, '_', sat_ID,'_',date_curr))
% 
% reci = Xsp3_eci(1:3,:); veci = Xsp3_eci(4:6,:); yyyy = yyyy_calc; mon = mon_calc; day_data = day_calc;
% hh = hh_calc; mm = mm_calc; ss = ss_calc; recef = data_.sp3_p;
% time_vec_uni = time_uni;
% 
% jdgps_sec_data_prop = 86400*GREGORIANtoJD_vector(yyyy,mon,day_data) + 3600*hh + 60*mm + ss; %%%%%%%%%%%%
% jdutc_sec_data_prop = jdgps_sec_data_prop - leap_sec + 19;
% 
% 
% X_true_aug = interp1(ode_step + parameters.jd_init, X1_prop,jdutc_sec_data_prop,'spline');
% X_true_samp = X_true_aug';
% % [t,X1] = ode45(@(t,x)propagator_truth(t, x, parameters), time_pred_utc, X_true_aug(:,1), options);
% % X_true_aug = X1';
% %
% % X_true_samp = interp1(time_pred_utc, X_true_aug(1:6,:)', sod_sec_sp3 - sod_sec_sp3(1), 'spline');
% % X_true_samp = X_true_samp';
% 
% err_x_drag = X_true_samp(1:6,:) - [reci;veci];
% 
% parameters_temp = parameters;
% parameters_temp.sun       = 0;
% parameters_temp.moon      = 0;
% parameters_temp.tides     = 0;
% parameters_temp.relativity= 0;
% flag_tides.SolidEarthTides = 0;
% flag_tides.OceanTides = 0;
% deg_grav = 0;
% parameters_temp.deg_grav = deg_grav;
% parameters_temp.ord_grav = deg_grav;
% parameters_temp.Cnm = Cbar(1:deg_grav+1, 1:deg_grav+1);
% parameters_temp.Snm = Sbar(1:deg_grav+1, 1:deg_grav+1);
% %             parameters.vec_est = [1:9];
% time_prop_utc_ekf = jdutc_sec_data - jdutc_sec_data(1);
% for ii = 1:numel(time_prop_utc_ekf)
%     X_full(vec_est,ii) = Xs_est(:,ii);
%     [~,Cd_prop(ii),rho_prop(ii),rho_nom_prop(ii),a_total(ii), a_drag(:,ii), a_srp(:,ii),a_earth_rad(:,ii), theta_calc(:,ii), phi_calc(ii),rot_ECI2ECEF] = ...
%         propagator_dmc_real(time_prop_utc_ekf(ii),X_full(:,ii),parameters_temp);
% end