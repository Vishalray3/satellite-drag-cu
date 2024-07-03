%% function
function [Cd_total,frac,P_o] = cd_body_exp(angle_step,parameters,trig)
E = angle_step;
e = parameters.e;
a_sma = parameters.a_sma;
inc = parameters.inc;
raan = parameters.raan;
w_arg = parameters.w_arg;
mu_e = parameters.mu_e;
str_special = parameters.str_special;
u_arg = parameters.u_arg;
R = parameters.R;
Kl = parameters.Kl;
M_s = parameters.M_s;
Re = parameters.Re;
Fwind = parameters.Fwind;
Ar = parameters.Ar;
A = parameters.Area_plates;
Tw = parameters.Tw;
Alpha = parameters.Alpha;
N = parameters.Order_b;
Ai = parameters.area_vec;
flag_axis = parameters.flag_axis;
shape_model = parameters.shape_model;
amu = parameters.amu;
oxa = parameters.atm_mass(2);
k_b = parameters.k_b;
phi = parameters.phi;
theta = parameters.theta;

sin_ano = sqrt(1-e^2)*sin(E)./(1-e*cos(E));
cos_ano = (cos(E) - e)./(1-e*cos(E));
true_ano = atan2d(sin_ano,cos_ano);
u_arg = true_ano;
M_ano = E-e*sin(E);
t = M_ano/parameters.n_mean;
time_et = parameters.time_et + t;
[r_eci, v_eci] = coe2rv(a_sma,e,inc,raan,w_arg,true_ano,mu_e,str_special, u_arg);
X_state = [r_eci;v_eci];

[~, ~, M, rho_o, T] = density_output(t, parameters.epoch, parameters.doy, time_et, X_state, parameters.sun_pos, parameters.eps, parameters);

Vi = Fwind*norm(v_eci);

Ri = R/M;
S = Vi./sqrt(2.*Ri.*T);

P_o = rho_o/(oxa*amu)*k_b*T;
m_r = M./M_s;

frac = Kl.*P_o./(1+Kl.*P_o);
if strcmp(shape_model, 'plate_quasi')
    Cd_total = single_plate_quasi_cd(phi, theta, S, Ai, m_r, Ri, Tw, Vi, A, Ar, Alpha, frac, flag_axis);
elseif strcmp(shape_model, 'plate_dria')
    Cd_total = single_plate_dria_cd(phi, theta, S, Ai, m_r, Ri, Tw, Vi, A, Ar, Alpha, frac, flag_axis);
elseif strcmp(shape_model,'sphere') || strcmp(shape_model,'sphere_jac')
    Cd_total = sphere_dria_cd(S, m_r(1), Ri, Tw, Vi, Alpha,frac,Kl,P_o,shape_model);
end
m = [0:N]';
if trig ==1
    Cd_total = repmat(Cd_total, numel(m),1).*repmat(cos(m*theta),1,numel(Cd_total));% Cd_total = Cd_total.*cos(m*theta).*cos(n*phi);
elseif trig ==2
    Cd_total = repmat(Cd_total, numel(m),1).*repmat(sin(m*theta),1,numel(Cd_total)); %Cd_total = Cd_total.*sin(m*theta).*cos(n*phi);
    % elseif trig == 3
    %     Cd_total = Cd_total.*cos(m*theta).*sin(n*phi);
    % elseif trig == 4
    %     Cd_total = Cd_total.*sin(m*theta).*sin(n*phi);
end
end