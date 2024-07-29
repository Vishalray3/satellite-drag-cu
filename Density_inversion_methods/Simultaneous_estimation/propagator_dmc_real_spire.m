%% propagator
function [Xdot,Cd_est,rho, rho_nom, a_total, a_drag, a_srp_ntw,a_earthrad_ntw,theta_in, phi_in, rot_ECI2ECEF, a_grav_eci, a_sun_eci, a_moon_eci,...
    a_srp_eci, a_earthrad_eci, t_ind] = propagator_dmc_real_spire(t,X,parameters)
% %#codegen
% coder.extrinsic('matlab')
%% Parameters
N_st = parameters.N_st;
N_nom = parameters.N_nom;
Order_o = parameters.Order_o;
Order_b = parameters.Order_b;
Af_ind = parameters.Af_ind;
Bf_ind = parameters.Bf_ind;
Cf_ind = parameters.Cf_ind;
Df_ind = parameters.Df_ind;
t_et  = parameters.time_et;
tym = parameters.tym;
omega_e = parameters.omega_e;
epoch = parameters.epoch;                    %%%%%%% REMOVE THIS
doy = parameters.doy;                   %%%%%%% REMOVE THIS
eps = parameters.eps;
flag_srp = parameters.flag_srp;
estimated_coeff = parameters.estimated_coeff;
vec_est = parameters.vec_est;
omega = parameters.omega;
zeta = parameters.zeta;
tau_inv = parameters.tau_inv;
tau_inv_cd = parameters.tau_inv_cd;
c1 = parameters.c1;
c2 = parameters.c2;
c3 = parameters.c3;
c_cd = parameters.c_cd;
Cd_est = parameters.Cd_est;
Cnm = parameters.Cnm;
Snm = parameters.Snm;
mu_e = parameters.mu_e;
% theta_max = parameters.theta_max;
q_bo = parameters.q_bo;
Cr_est = parameters.Cr_est;
eop = parameters.eop;
NUT1980Info = parameters.nut80;
leap_sec = parameters.leapSec;
c_cr = parameters.c_cr;
% tauinv_accn = parameters.tauinv_accn;
% tauinv_acct = parameters.tauinv_acct;
% tauinv_accw = parameters.tauinv_accw;
% c_accn = parameters.c_accn;
% c_acct = parameters.c_acct;
% c_accw = parameters.c_accw;
% scale_acc = parameters.scale_acc;
Cr_erp = parameters.Cr_erp;
% rollmat = parameters.roll;
% pitchmat = parameters.pitch;
% yawmat = parameters.yaw;
jd_init = parameters.jd_init;

X_state = X(1:6);

% theta_in = zeros(numel(theta(:,1)),1);
accel = zeros(3,1);
acc_new = zeros(3,1);
Fpos = zeros(3,3);
Fvel = zeros(3,3);
Fpos_new = zeros(3,3);
Fvel_new = zeros(3,3);
coder.varsize('Fa_par')
coder.varsize('Fpar_new')
dut1 = 0;
xp = 0;
yp = 0;
rot_ECI2ECEF = zeros(3,3);
%% Interpolation


[~,t_ind] = min(abs(tym - t));
% if strcmp(flag_drag,'Bod')
%     %           phi_in = Input(2,t_ind);
%     %           theta_in = Input(3,t_ind);
%     mu_e = Force(1).Model.mu_e;
%     coe = rv2coe_E(X_state(1:3),X_state(4:6),mu_e);
%     theta_in = coe(7); % for circulr orbits, take the argument of latitude
%     phi_in = coe(7);
% end
% if strcmp(flag_drag,'Bod_orb')
%
%     mu_e = Force(1).Model.mu_e;
%     coe = rv2coe_E(X_state(1:3),X_state(4:6),mu_e);
%     theta_in = coe(7); % for circulr orbits, take the argument of latitude
%     phi_in = coe(7);
% end

if t_ind ~= numel(tym)
    if tym(t_ind) ~= tym(t_ind+1)
        t2 = tym(t_ind+1);
        
    else
        t2 = tym(t_ind+2);
    end
    t1 = tym(t_ind);
    time_et = (t-t2)*(t_et(t_ind+1) - t_et(t_ind))/(t2-t1) + t_et(t_ind+1);
else
    t2 = tym(t_ind);
    if tym(t_ind) ~= tym(t_ind-1)
        t1 = tym(t_ind-1);
        
    else
        t1 = tym(t_ind-2);
    end
    time_et = (t-t2)*(t_et(t_ind) - t_et(t_ind-1))/(t2-t1) + t_et(t_ind);
end

time_jd = jd_init + t;
q_inst = q_bo(:,t_ind);


X_f = X(N_nom+1:N_nom+parameters.N_f,1);
% if strcmp(flag_drag,'Orb') == 1
%     %%%%%%%%%%%%%%%%% Orbit Coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     theta_in = 0;
%     mu_e = Force(1).Model.mu_e;
%     coe = rv2coe_E(X_state(1:3),X_state(4:6),mu_e);
%     phi_in = coe(7);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Drag coefficient calculation

% %%%%%%
% if theta_in < 0
%     theta_in = -theta_in;
% end
% theta_in = theta_in/theta_max;
% %%%%%%



%% ECEF to ECI rotation matrix

%[rot_ECI2ECEF,dut1,xp,yp] = time2rotmat(eop, time_var,X_state, parameters.eqeterms, parameters.nut80, [0;0;0]); % parameters.rot_ECI2ECEF(:,:,t_ind);
[~,~,rot_ECI2ECEF,dut1,xp,yp] = PosVelConvert(X_state(1:3)'/1000, X_state(4:6)'/1000,time_jd,...
    'J2K2ECF','106terms', leap_sec, NUT1980Info, eop);
% rot_ECI2ECEF = cspice_pxform( 'J2000', 'ITRF93', time_var );
parameters.TEI = rot_ECI2ECEF;

%% ECI TO NTW rotation matrix
rhat = X_state(1:3)/norm(X_state(1:3));
vhat = X_state(4:6)/norm(X_state(4:6));
That = vhat;
What = cross(vhat,rhat);
What = What/norm(What);
Nhat = cross(That,What);
Nhat = Nhat/norm(Nhat);
ECI2NTW = [ Nhat'; That'; What'];

%% ECI TO BODY
rot_ECItoVVLH(3,:) = -X_state(1:3)'/norm(X_state(1:3));
yhat               = cross(X_state(4:6)/norm(X_state(4:6)),X_state(1:3)/norm(X_state(1:3)));
rot_ECItoVVLH(2,:) = yhat'/norm(yhat);
rot_ECItoVVLH(1,:) = cross(rot_ECItoVVLH(2,:),rot_ECItoVVLH(3,:));
rot_ECItoVVLH(1,:) = rot_ECItoVVLH(1,:)/norm(rot_ECItoVVLH(1,:));

rot_VVLHtoSBF = quat2rotmat(q_inst);

SBF2ECI = rot_ECItoVVLH'*rot_VVLHtoSBF';
SBF2ECEF = rot_ECI2ECEF*SBF2ECI;

% sp_angle = parameters.sada_angle(t_ind);
% parameters.area_vec(:,7) = [0; cosd(sp_angle); sind(sp_angle)];
% parameters.area_vec(:,8) = -[0; cosd(sp_angle); sind(sp_angle)];
n_hat = SBF2ECI*parameters.area_vec;           % body to eci of plate vectors

%% Calculate theta and phi
omega_vec = [0;0;1]*omega_e;
corotwind_eci = cross(omega_vec, X_state(1:3));
vrel_eci = X_state(4:6) - corotwind_eci;
vrel_sbf = SBF2ECI'*vrel_eci;
theta_in = atan2d(vrel_sbf(2),vrel_sbf(1));         % angle in x-y plane
phi_in = atan2d(vrel_sbf(3),sqrt(vrel_sbf(1)^2+vrel_sbf(2)^2));  % angle with z-axis
%% Force models
N_est = 6; %numel(vec_est);
Ncr = 0;
Ba = zeros(6,1);
Fpar_dot = [];
Xpar_dot = [];


if parameters.tides
    tai_sec = time_jd + leap_sec;
    tt_sec = tai_sec + 32.184;
    jdTT_ref = (tt_sec/86400 - 2451545)/36525;
    jd_ut1 = (time_jd + dut1)/86400 - 2451545;
    [Cnm, Snm] = tides_fes(parameters.Cnm, parameters.Snm, time_et, jdTT_ref, jd_ut1, rot_ECI2ECEF, parameters.moon_pos(:,t_ind),...
        parameters.sun_pos(:,t_ind), parameters.Re, mu_e,parameters.mu_moon, parameters.mu_sun, xp, yp, parameters.flag_tides, ...
        parameters.coeff0, parameters.coeff1, parameters.coeff2, parameters.doo_coeff, parameters.Cnmp_ot, parameters.Snmp_ot, ...
        parameters.Cnmm_ot, parameters.Snmm_ot, parameters.deg_do, parameters.ord_do, parameters.pole_mean);
end


Fmod1          = HighGrav1(parameters.deg_grav, parameters.ord_grav, mu_e, parameters.Re, Cnm, Snm, rot_ECI2ECEF);
[accel,Fpos,Fvel,Fa_par]  = Fmod1.compAccelPar(X_state,t);

a_grav_eci = accel;
if parameters.sun
    Fmod2      = ThreeBody(parameters.mu_sun, parameters.sun_pos(:,t_ind));
    [acc_new,Fpos_new,Fvel_new,~]  = Fmod2.compAccelPar(X_state,t);
    accel = accel + acc_new;
    Fpos = Fpos + Fpos_new;
    Fvel = Fvel + Fvel_new;
end

a_sun_eci = acc_new;
if parameters.moon
    Fmod3      = ThreeBody(parameters.mu_moon, parameters.moon_pos(:,t_ind));
    [acc_new,Fpos_new,Fvel_new,~]  = Fmod3.compAccelPar(X_state,t);
    accel = accel + acc_new;
    Fpos = Fpos + Fpos_new;
    Fvel = Fvel + Fvel_new;
end
a_moon_eci = acc_new;

if parameters.srp
    SRP_am    = parameters.srp_area/parameters.mass;
    if strcmp(flag_srp, 'Cball')
        if estimated_coeff.Cr
            Ncr = 1;
            Cr_est = X(7,1);
            
        end
        Fmod41      = SRP(SRP_am, parameters.AU, Cr_est, parameters.P_sun, parameters.sun_pos(:,t_ind), parameters.Re, parameters.R_sun);
        [acc_new,Fpos_new,Fvel_new,Fpar_new]  = Fmod41.compAccelPar(X_state,t);
        
    elseif strcmp(flag_srp, 'Three')
        if estimated_coeff.Cr
            Ncr = 3;
            A0 = X(7,1);
            A1 = X(8,1);
            A2 = X(9,1);
            c_cr = [c_cr;c_cr;c_cr];
        end
        Fmod42     = SRP_3constants(SRP_am, parameters.AU, Cr_est, parameters.P_sun, parameters.sun_pos(:,t_ind),...
            parameters.Re, parameters.R_sun, parameters.earth_vel(:,t_ind),A0,A1,A2);
        [acc_new,Fpos_new,Fvel_new,Fpar_new]  = Fmod42.compAccelPar(X_state,t);
        
    elseif strcmp(flag_srp, 'Panel')
        %         ang_pan = parameters.angle_panels(t_ind,:,:);
        %         zero_mat = zeros(1,1,numel(ang_pan));
        %         one_mat = ones(1,1,numel(ang_pan));
        %         rot_pan2bod = [cosd(ang_pan), zero_mat, -sind(ang_pan); zero_mat,one_mat,zero_mat; sind(ang_pan), zero_mat, cosd(ang_pan)];
        %         area_vec = reshape(parameters.area_vec,3,1,9);
        %         area_vec_b = squeeze(pagemtimes(rot_pan2bod,area_vec));
        %         SBF2ECI = parameters.rot_SBF2ECI(:,:,t_ind);
        %         n_hat = SBF2ECI*area_vec_b;
        %         SBF2ECI = parameters.rot_SBF2ECI(:,:,t_ind);
        if estimated_coeff.Cr
            Ncr = 1;
            Cr_est = X(7,1);
        end
        Fmod43      = SRP_GPM(parameters.mass, parameters.AU, parameters.rho_spec, parameters.rho_diff,parameters.P_sun,...
            parameters.sun_pos(:,t_ind), parameters.Re, parameters.R_sun , n_hat,parameters.Area_plates,numel(parameters.rho_spec),Cr_est);
        [acc_new,Fpos_new,Fvel_new,Fpar_new]  = Fmod43.compAccelPar(X_state,t);
    end
    accel = accel + acc_new;
    Fpos = Fpos + Fpos_new;
    Fvel = Fvel + Fvel_new;
    a_srp_ntw = ECI2NTW*acc_new;
    a_srp_eci = acc_new;
    
    if estimated_coeff.Cr
        Fa_par = [Fa_par, Fpar_new];
        Fpar_dot = [Fpar_dot; zeros(Ncr,N_st)];
        Xpar_dot = [Xpar_dot; zeros(Ncr,1)];
        N_est = N_est+Ncr;
        %         Ba = [Ba;zeros(Ncr,1)];
        Ba = [Ba;c_cr];
    end
    
end


if parameters.drag
    %% density calculation
    
    if estimated_coeff.Cd
        Cd_est = Cd_est+ X(N_nom,1);
        
        n_grid = repmat([0:Order_o],Order_b+1,1);
        m_grid = repmat([0:Order_b]',1,Order_o+1);
        Af_arg = cosd(n_grid*phi_in).*cosd(m_grid*theta_in(1));
        Bf_arg = cosd(n_grid*phi_in).*sind(m_grid*theta_in(1));
        Cf_arg = sind(n_grid*phi_in).*cosd(m_grid*theta_in(1));
        Df_arg = sind(n_grid*phi_in).*sind(m_grid*theta_in(1));
        % nominal entries corresponing to the indices
        Af_filter = Af_arg(Af_ind); Bf_filter = Bf_arg(Bf_ind); Cf_filter = Cf_arg(Cf_ind); Df_filter = Df_arg(Df_ind);
        % vectors from those entries
        argf_vec = [Af_filter(:); Bf_filter(:); Cf_filter(:); Df_filter(:)];
        argf_vec = argf_vec(2:end);
        
        Cd_est = Cd_est + sum(X_f.*argf_vec);
        [rho_nom, delrho] = density_output(time_jd, X_state, parameters.sun_pos(:,t_ind), eps, parameters);
        bff_coeff = [zeros(Order_b+1,1),zeros(Order_b+1,1)];
    elseif estimated_coeff.CdTrue
        parameters.t = t;
        parameters.time_et = time_jd;
        parameters.theta = pi/180*theta_in';
        parameters.X_state = X_state;
        parameters.phi = pi/180*phi_in';
        parameters.v_sbf = [];%parameters.v_sbf_mat(t_ind,:)';
        [rho_nom,delrho,Cd_est,Aref,Cl_est, Cd_ads, Cd_s] = cd_truth_bff(time_jd,parameters);
    elseif estimated_coeff.CdDiscrete
        [rho_nom, delrho] = density_output(time_jd, X_state, parameters.sun_pos(:,t_ind), eps, parameters);
        Cd_est = parameters.Cd_discrete(t_ind);
        Cd_ads = 0;
        Cd_s = 0;
    elseif estimated_coeff.CdDiscrete_est
        [rho_nom, delrho] = density_output(time_jd, X_state, parameters.sun_pos(:,t_ind), eps, parameters);
        Cd_est = parameters.Cd_discrete(t_ind)*X(N_nom);
    else
        [rho_nom, delrho] = density_output(time_jd, X_state, parameters.sun_pos(:,t_ind), eps, parameters);
        bff_coeff = [zeros(Order_b+1,1),zeros(Order_b+1,1)];
    end
    if estimated_coeff.rho_DMC
        rho = rho_nom*(1 + X(6+Ncr+1)+X(6+Ncr+3));
    else
        rho = rho_nom;
    end
    Fmod5      = Drag_atm(parameters.mass, parameters.drag_area, Cd_est, omega_vec, rho, delrho, [0;0;0]);
    [acc_new,Fpos_new,Fvel_new,Fpar_new]  = Fmod5.compAccelPar(X_state,t);
    a_drag = acc_new;
    accel = accel + acc_new;
    Fpos = Fpos + Fpos_new;
    Fvel = Fvel + Fvel_new;
    if parameters.lift
        Fmodl      = Lift_atm(parameters.mass, parameters.drag_area, Cl_est, omega_vec, rho, delrho, [0;0;0],n_hat,numel(parameters.Area_plates));
        [acc_new1,Fpos_new,Fvel_new,Fpar_new]  = Fmodl.compAccelPar(X_state,t);
        a_lift = norm(acc_new1);
        accel = accel + acc_new1;
        Fpos = Fpos + Fpos_new;
        Fvel = Fvel + Fvel_new;
        acc_new = acc_new + acc_new1;
    else
        a_lift = [0;0;0];
    end
    if estimated_coeff.rho_DMC
        Fdmc_rho = [acc_new/rho*rho_nom, zeros(3,1), acc_new/rho*rho_nom];
        Fa_par = [Fa_par, Fdmc_rho];
        Flast_dmc_rho = [zeros(1,N_est), 0, 1,0, zeros(1,N_st-N_est-3); zeros(1,N_est), -omega^2, -2*zeta*omega,0, zeros(1,N_st-N_est-3);...
            zeros(1,N_est),0,0,-tau_inv,zeros(1,N_st-N_est-3)];
        Fpar_dot = [Fpar_dot; Flast_dmc_rho];
        Ba = [Ba;c1;c2;c3];
        Xpar_dot = [Xpar_dot; 0*X(6+Ncr+1)+X(6+Ncr+2); -omega^2*X(6+Ncr+1)-2*zeta*omega*X(6+Ncr+2);-tau_inv*X(6+Ncr+3)];
        N_est = N_est + 3;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if estimated_coeff.Cd
        Fa_par = [Fa_par, Fpar_new];
        if numel(argf_vec) > 0
            %         Fcd2 = Fcd.*argf_vec(vec_est((N_nom+1):end)-N_nom)';
            Fpar_new = Fpar_new.*argf_vec';
            Fa_par = [Fa_par, Fpar_new];
        end
        Fpar_dot = [Fpar_dot;zeros(1,N_nom-1),-tau_inv_cd, zeros(1,N_st-N_nom); zeros(N_st-N_nom,N_st)];
        Ba = [Ba;c_cd;zeros(N_st-N_nom,1)];
        Xpar_dot = [Xpar_dot;-tau_inv_cd*X(N_nom);zeros(N_st-N_nom,1)];
        N_est = N_est + (N_st-N_nom+1);
    elseif estimated_coeff.CdDiscrete_est
        Fpar_new = Fpar_new*parameters.Cd_discrete(t_ind);
        Fa_par = [Fa_par, Fpar_new];
        Fpar_dot = [Fpar_dot;zeros(1,N_nom-1),-tau_inv_cd, zeros(1,N_st-N_nom); zeros(N_st-N_nom,N_st)];
        Ba = [Ba;c_cd;zeros(N_st-N_nom,1)];
        Xpar_dot = [Xpar_dot;-tau_inv_cd*X(N_nom);zeros(N_st-N_nom,1)];
        N_est = N_est + (N_st-N_nom+1);
    end
end

if parameters.relativity
    Fmod6      = Relativity(parameters.c,parameters.mu_e);
    [acc_new,Fpos_new,Fvel_new,~]  = Fmod6.compAccelPar(X_state,t);
    accel = accel + acc_new;
    Fpos = Fpos + Fpos_new;
    Fvel = Fvel + Fvel_new;
end

if parameters.earth_rad
    if estimated_coeff.Cr_erp
        Cr_erp = X(N_est+1);
        N_est = N_est+1;
    end
    if strcmp(parameters.flag_earthrad, 'Knocke')
    sod = mod(epoch+t, 86400);
    Fmod7 = bwalbedo(parameters.year - 2000,doy,sod,parameters.sun_pos(:,t_ind),parameters.area_vec,parameters.Area_plates,parameters.mass,...
        parameters.rho_spec,parameters.rho_diff,parameters.rho_spec_ir,parameters.rho_diff_ir,parameters.TEI',SBF2ECEF, Cr_erp);
    else
       Fmod7 = ceres_erp_fast(parameters.sun_pos(:,t_ind),n_hat,parameters.Area_plates,parameters.mass,parameters.rho_spec,parameters.rho_diff,...
            parameters.rho_spec_ir,parameters.rho_diff_ir,parameters.TEI', parameters.Psw, parameters.Plw, parameters.lon, parameters.lat, ...
            Cr_erp, parameters.Re); 
    end
    [acc_new,Fpos_new,Fvel_new,Fpar_new]  = Fmod7.compAccelPar(X_state,t);
    accel = accel + acc_new;
    Fpos = Fpos + Fpos_new;
    Fvel = Fvel + Fvel_new;
    a_earthrad_ntw = ECI2NTW*acc_new;
    a_earthrad_eci = acc_new;
    if estimated_coeff.Cr_erp
        Fa_par = [Fa_par, Fpar_new];
        Fpar_dot = [Fpar_dot;zeros(1,N_st)];
        Xpar_dot = [Xpar_dot;0];
        Ba = [Ba;0];
    end
else
    a_earthrad_ntw = [0;0;0];
    a_earthrad_eci = acc_new;
end

if parameters.empirical
    %     ECI2NTW = eye(3);
    acc_new = ECI2NTW'*X(N_est+1:N_est+3)*scale_acc;
    accel = accel + acc_new;
    Fpard = diag([-tauinv_accn, -tauinv_acct, -tauinv_accw]);
    Fpard = ECI2NTW'*Fpard*ECI2NTW;
    Fpar_dot = [Fpar_dot; zeros(3, N_st-3), Fpard ];
    Fpar_new = ECI2NTW'*ECI2NTW*scale_acc;
    Fa_par = [Fa_par, Fpar_new];
    Fpar
    Xpar_dot = [Xpar_dot;-tauinv_accn*X(N_est+1);-tauinv_acct*X(N_est+2);-tauinv_accw*X(N_est+3)];
    Ba = [Ba; c_accn; c_acct; c_accw];
end
Ba = Ba(vec_est);
a_total = norm(accel);
%% Jacobian calc
F = [zeros(3,3), eye(3),zeros(3,N_st-6);
    Fpos, Fvel,Fa_par;
    Fpar_dot];

%% Stm dynamics

if strcmp(parameters.case_est,'EKF')
    F = F(vec_est,vec_est);
    N_ac = numel(vec_est);
    stm = reshape(X(N_st+1:N_st+N_ac^2),N_ac,N_ac);
    stm_dot = F*stm;
    X_state_dot = [X_state(4:6);accel;Xpar_dot];
    stm_q_dot = stm*(Ba*Ba')*stm';
    %% Integrating the smoothing vector b
    b = X(N_st+2*N_ac^2+1:end);
    % X_state_dot(vec_est) - F*X(vec_est)
    bdot = F*b + X_state_dot(vec_est) - F*X(vec_est);
    %% Forming the dynamics vector
    Xdot = [X_state_dot;stm_dot(:);stm_q_dot(:);bdot];                            % add dynamics of other states if there
elseif strcmp(parameters.case_est,'EDR')
    X_state_dot = [X_state(4:6);accel;Xpar_dot];
    Xdot = X_state_dot;
else
    stm = reshape(X(N_st+1:N_st+N_st^2),N_st,N_st);
    stm_dot = F*stm;
    X_state_dot = [X_state(4:6);accel;Xpar_dot];
    Xdot = [X_state_dot;stm_dot(:)];
end
end


