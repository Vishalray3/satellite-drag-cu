%% propagator
function [Xdot,Cd_est,rho,theta_in,phi_in,bff_coeff,S,m_r,frac,r_ads] = propagator_dmc_log(t,X,parameters)
%% Parameters
N_st = parameters.N_st;
N_nom = parameters.N_nom;
Order_o = parameters.Order_o;
Order_b = parameters.Order_b;
Af_ind = parameters.Af_ind;
Bf_ind = parameters.Bf_ind;
Cf_ind = parameters.Cf_ind;
Df_ind = parameters.Df_ind;
time_prop = parameters.time_prop;
tym = parameters.tym;
omega_e = parameters.omega_e;
epoch = parameters.epoch;
doy = parameters.doy;
eps = parameters.eps;
Force = parameters.Force;
flag_srp = parameters.flag_srp;
flag_drag = parameters.flag_drag;
sun_positions = parameters.sun_pos;
moon_positions = parameters.moon_pos;
earth_velocity = parameters.earth_vel;
bo_mod = parameters.bo_mod;
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
T_orb = parameters.T_orb;
% theta_max = parameters.theta_max;
theta = parameters.theta;
Cr_est = parameters.Cr_est;
if ismember('Cr',estimated_coeff)
    if strcmp(flag_srp, 'Cball')
        Cr_est = X(7,1);
    elseif strcmp(flag_srp,'Three')
        A0_est = X(7,1);
        A1_est = X(8,1);
        A2_est = X(9,1);
    end
end

accel = zeros(3,1);
X_state = X(1:6);
Fpos = zeros(3,3);
Fvel = zeros(3,3);

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
    time_var = (t-t2)*(time_prop(t_ind+1) - time_prop(t_ind))/(t2-t1) + time_prop(t_ind+1);
    ind_theta = (theta(:,t_ind+1)./abs(theta(:,t_ind+1))).*(theta(:,t_ind)./abs(theta(:,t_ind))) > 0;
    theta_in(ind_theta) = (t-t2)*(theta(ind_theta,t_ind+1) - theta(ind_theta,t_ind))/(t2-t1) + theta(ind_theta,t_ind+1);
    theta_in(~ind_theta) = theta(~ind_theta,t_ind);
    mu_e = Force(1).Model.mu_e;
    coe = rv2coe_E(X_state(1:3),X_state(4:6),mu_e);
    phi_in = coe(7); 
else
    t2 = tym(t_ind);
    if tym(t_ind) ~= tym(t_ind-1)
        t1 = tym(t_ind-1);
        
    else
        t1 = tym(t_ind-2);
    end
    time_var = (t-t2)*(time_prop(t_ind) - time_prop(t_ind-1))/(t2-t1) + time_prop(t_ind);
    ind_theta = (theta(:,t_ind)./abs(theta(:,t_ind))).*(theta(:,t_ind-1)./abs(theta(:,t_ind-1))) > 0;
    theta_in(ind_theta) = (t-t2)*(theta(ind_theta,t_ind) - theta(ind_theta,t_ind-1))/(t2-t1) + theta(ind_theta,t_ind);
    theta_in(~ind_theta) = theta(~ind_theta,t_ind);
    mu_e = Force(1).Model.mu_e;
    coe = rv2coe_E(X_state(1:3),X_state(4:6),mu_e);
    phi_in = coe(7);
end
X_f = X(N_nom+1:N_st,1);
% if strcmp(flag_drag,'Orb') == 1
%     %%%%%%%%%%%%%%%%% Orbit Coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     theta_in = 0;
%     mu_e = Force(1).Model.mu_e;
%     coe = rv2coe_E(X_state(1:3),X_state(4:6),mu_e);
%     phi_in = coe(7);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sun_pos = sun_positions(:,t_ind);
moon_pos = moon_positions(:,t_ind);
earth_vel = earth_velocity(:,t_ind);

%% SRP calculations
ang_pan = parameters.angle_panels(t_ind,:,:);
zero_mat = zeros(1,1,numel(ang_pan));
one_mat = ones(1,1,numel(ang_pan));
rot_pan2bod = [cosd(ang_pan), zero_mat, -sind(ang_pan); zero_mat,one_mat,zero_mat; sind(ang_pan), zero_mat, cosd(ang_pan)];
area_vec = reshape(parameters.area_vec,3,1,9);
area_vec_b = squeeze(pagemtimes(rot_pan2bod,area_vec));
SBF2ECI = parameters.rot_SBF2ECI(:,:,t_ind);
n_hat = SBF2ECI*area_vec_b;                          % body to eci of plate vectors
% SBF2ECI = parameters.rot_SBF2ECI(:,:,t_ind);
% n_hat = SBF2ECI*parameters.area_vec;           % body to eci of plate vectors
%% Drag coefficient calculation

% %%%%%%
% if theta_in < 0
%     theta_in = -theta_in;
% end
% theta_in = theta_in/theta_max;
% %%%%%%

if ismember('Cd',estimated_coeff)
    Cd_est = Cd_est+ X(N_nom,1);
    
    n_grid = repmat([0:Order_o],Order_b+1,1);
    m_grid = repmat([0:Order_b]',1,Order_o+1);
    Af_arg = cosd(n_grid*phi_in).*cosd(m_grid*theta_in);
    Bf_arg = cosd(n_grid*phi_in).*sind(m_grid*theta_in);
    Cf_arg = sind(n_grid*phi_in).*cosd(m_grid*theta_in);
    Df_arg = sind(n_grid*phi_in).*sind(m_grid*theta_in);
    % nominal entries corresponing to the indices
    Af_filter = Af_arg(Af_ind); Bf_filter = Bf_arg(Bf_ind); Cf_filter = Cf_arg(Cf_ind); Df_filter = Df_arg(Df_ind);
    % vectors from those entries
    argf_vec = [Af_filter(:); Bf_filter(:); Cf_filter(:); Df_filter(:)];
    argf_vec = argf_vec(2:end);
    
    Cd_est = Cd_est + sum(X_f.*argf_vec);
end

%% ECEF to ECI rotation matrix
% rot_ECI2ECEF = time2rotmat(eop, time_var,X_state);
rot_ECI2ECEF = cspice_pxform( 'J2000', 'ITRF93', time_var );
parameters.TEI = rot_ECI2ECEF;
%% density calculation
if ismember('CdTrue',estimated_coeff)
    parameters.t = t;
    parameters.time_et = time_var;
    parameters.theta = pi/180*theta_in;
    parameters.X_state = X_state;
    parameters.phi = pi/180*parameters.phi(:,t_ind)';
    [rho_nom,delrho,Cd_est, bff_coeff,S,m_r,frac,r_ads] = cd_truth_bff(t,parameters);
else
    [rho_nom, delrho] = density_output(t, epoch, doy, time_var, X_state, sun_pos, eps, parameters);
    bff_coeff = [zeros(Order_b+1,1),zeros(Order_b+1,1)];
    
end
omega_vec = [0;0;1]*omega_e;
%% Force models
if ismember('rho_DMC',estimated_coeff) && ~ismember('Cr',estimated_coeff)
    rho = rho_nom*(1+X(7)+X(9));
elseif ismember('rho_DMC',estimated_coeff) && ismember('Cr',estimated_coeff)
    rho = rho_nom*(1+X(8)+X(10));
else
    rho = rho_nom;
end

Mjd_UTC = parameters.Mjd_UTC0 + t/86400;
index_dut = parameters.index_dut + floor(t/86400);
UT1_UTC = parameters.dut(index_dut) + (parameters.dut(index_dut+1)-parameters.dut(index_dut))*mod(t/86400,1) ;
x_pole = parameters.x_pole(index_dut) + (parameters.x_pole(index_dut+1)-parameters.x_pole(index_dut))*mod(t/86400,1) ;
y_pole = parameters.y_pole(index_dut) + (parameters.y_pole(index_dut+1)-parameters.y_pole(index_dut))*mod(t/86400,1) ;
[Cnm, Snm] = tides(parameters.Cnm, parameters.Snm, time_var, Mjd_UTC, UT1_UTC, rot_ECI2ECEF, moon_pos, sun_pos, parameters.Re, mu_e,...
                   parameters.mu_m, parameters.mu_s, x_pole, y_pole, parameters.flag_tides);
Force(1).Model.Cbar  = Cnm;
Force(1).Model.Sbar  = Snm;
% Force(1).Model.et = time_et;
Force(1).Model.TEI = rot_ECI2ECEF;

Force(2).Model.Cd        = Cd_est;
Force(2).Model.rho       = rho;
Force(2).Model.delrho    = delrho;
Force(2).Model.omega_vec = omega_vec;
Force(2).Model.X_f       = X_f;

Force(3).Model.X_sun = sun_pos;
if strcmp(flag_srp, 'Cball')
    Force(3).Model.Cr    = Cr_est;
elseif strcmp(flag_srp, 'Three')
    Force(3).Model.A0    = A0_est;
    Force(3).Model.A1    = A1_est;
    Force(3).Model.A2    = A2_est;
    Force(3).Model.earth_vel = earth_vel;
elseif strcmp(flag_srp, 'Panel')
    Force(3).Model.Cr    = Cr_est;
    Force(3).Model.n_hat = n_hat;
end

Force(4).Model.X_third  = sun_pos;
Force(5).Model.X_third = moon_pos;
%% Dynamics calculation
N_acc = numel(Force);

for i=1:N_acc                                            % position and velocity dynamics
    [accel_new,Fpos_new,Fvel_new] = Force(i).Model.compAccelPar(X_state, t);   % acceleration
    accel = accel + accel_new;
    Fpos = Fpos + Fpos_new;
    Fvel = Fvel + Fvel_new;
    %     Fpos = Fpos + Force(i).Model.compParPos(X_state, t);    % position jacobian
    %     Fvel = Fvel + Force(i).Model.compParVel(X_state, t);               % velocity jacobian only drag
end

% Fcd = Force(4).Model.compParCd(X_state, t);
N_est = 6; %numel(vec_est);
B_st = [zeros(6,1)];
if ismember('Cr', estimated_coeff)
    Fcr = Force(3).Model.compParCr(X_state, t);
    F_last_cr = zeros(1,N_st);
    cr_xdot = 0;
    N_p = 1;
    N_est = N_est+1;
    B_cr = 0;
else
    Fcr = [];
    F_last_cr = [];
    B_cr = [];
    cr_xdot = [];
    N_p = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember('rho_DMC',estimated_coeff)
    Fdmc_rho = [Force(2).Model.compAccelPar(X_state,t)/rho*rho_nom, zeros(3,1), Force(2).Model.compAccelPar(X_state,t)/rho*rho_nom];
    Flast_dmc_rho = [zeros(1,N_est), -tau_inv, 1,0, zeros(1,N_st-N_est-3); zeros(1,N_est), -omega^2, -2*zeta*omega,0, zeros(1,N_st-N_est-3);...
        zeros(1,N_est),0,0,0,zeros(1,N_st-N_est-3)];
    B_rho = [c1;c2;c3];
    dmc_rhodot = [-tau_inv*X(7+N_p)+X(8+N_p); -omega^2*X(7+N_p)-2*zeta*omega*X(8+N_p);0];
    N_est = N_est + 3;
else
    Fdmc_rho = [];
    Flast_dmc_rho = [];
    B_rho = [];
    dmc_rhodot = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember('Cd',estimated_coeff)
    Fcd = Force(2).Model.compParCd(X_state, t);
    if numel(argf_vec) > 0
        %         Fcd2 = Fcd.*argf_vec(vec_est((N_nom+1):end)-N_nom)';
        Fcd2 = Fcd.*argf_vec';
    else
        Fcd2 = [];
    end
    F_last_cd = [zeros(1,N_nom-1),-tau_inv_cd, zeros(1,N_st-N_nom); zeros(N_st-N_nom,N_st)];
    B_cd = [c_cd;zeros(N_st-N_nom,1)];
    cd_xdot = [-tau_inv_cd*X(N_nom);zeros(N_st-N_nom,1)];
    N_est = N_est + (N_st-N_nom+1);
else
    Fcd = [];
    Fcd2 = [];
    F_last_cd = [];
    B_cd = [];
    cd_xdot = [];
end


Ba = [B_st;B_cr;B_rho;B_cd];
Ba = Ba(vec_est);
%% Jacobian calc
F = [zeros(3,3), eye(3),zeros(3,N_st-6);
    Fpos, Fvel,Fcr,Fdmc_rho, Fcd, Fcd2;
    F_last_cr; Flast_dmc_rho; F_last_cd];
F = F(vec_est,vec_est);
%% Stm dynamics
N_ac = numel(vec_est);
stm = reshape(X(N_st+1:N_st+N_ac^2),N_ac,N_ac);
stm_dot = F*stm;
stm_q_dot = stm*(Ba*Ba')*stm';
% if ismember('rho_DMC',estimated_coeff)
%     stm_q_dot = stm*(Ba*Ba')*stm';
% else
%     stm_q_dot = [];
% end

X_state_dot = [X_state(4:6);accel;cr_xdot;dmc_rhodot;cd_xdot];
%% Integrating the smoothing vector b
b = X(N_st+2*N_ac^2+1:end);
% X_state_dot(vec_est) - F*X(vec_est)
bdot = F*b + X_state_dot(vec_est) - F*X(vec_est);
%% Forming the dynamics vector
Xdot = [X_state_dot;stm_dot(:);stm_q_dot(:);bdot];                            % add dynamics of other states if there

end


