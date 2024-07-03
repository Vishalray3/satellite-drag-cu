%% function
function Cd_total = cd_bodf_exp_gpm(angle_step,parameters,trig)
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
order_body = parameters.Order_b;
order_orbit = parameters.Order_o;
Ai = parameters.area_vec;
flag_axis = parameters.flag_axis;
shape_case = parameters.shape_model;
amu = parameters.amu;
oxa = parameters.atm_mass(2);
k_b = parameters.k_b;

sin_ano = sqrt(1-e^2)*sin(E)./(1-e*cos(E));
cos_ano = (cos(E) - e)./(1-e*cos(E));
true_ano = atan2d(sin_ano,cos_ano);
u_arg = true_ano; % + w_arg
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

frac =  Kl.*P_o./(1+Kl.*P_o);
r_ads = sqrt(0.5*(1+Alpha*(4*Ri*Tw./Vi.^2 -1)));
r_ads = r_ads*ones(1,numel(A));
Alpha_s = 3.6*m_r./(1+m_r).^2;
r_s = sqrt(0.5*(1+Alpha_s*(4*Ri*Tw./Vi.^2 -1)));

if strcmp(flag_axis,'y')
    Ci = sqrt(Ai(1,:).^2+Ai(3,:).^2); 
    delta = atan2(Ai(3,:), Ai(1,:));
elseif strcmp(flag_axis, 'z')
    Ci = sqrt(Ai(1,:).^2+Ai(2,:).^2);
    delta = atan2(Ai(1,:), Ai(2,:));
end

An_ads = fourier_plate_dria(A, Ar, S, r_ads,Ci, order_body);
An_s = fourier_plate_dria(A, Ar, S, r_s, Ci,order_body);

An_plates = frac*An_ads + (1-frac)*An_s;
if strcmp(shape_case, 'plate_jac')
    An_plates = P_o/(1+Kl*P_o)^2*(An_ads-An_s);
end
order_vec_body = [0:order_body]';

cos_delta = cos(order_vec_body.*delta);
sin_delta = sin(order_vec_body.*delta);

An_total = An_plates.*cos_delta;
Bn_total = An_plates.*sin_delta;

ang_gpm = 3.75*pi/180;
cos_theta = cos(order_vec_body.*ang_gpm);
sin_theta = sin(order_vec_body.*ang_gpm);
% An_total(1,[1:5,8:9]) = sum(An_total(:,[1:5,8:9]).*cos_theta,1);
% Bn_total(1,[1:5,8:9]) = sum(Bn_total(:,[1:5,8:9]).*sin_theta,1);
An_total(1,[1:5]) = sum(An_total(:,[1:5]).*cos_theta,1);
Bn_total(1,[1:5]) = sum(Bn_total(:,[1:5]).*sin_theta,1);
An_total(2:end,[1:5]) = 0;
Bn_total(2:end,[1:5]) = 0;
An_total = sum(An_total,2);
Bn_total = sum(Bn_total,2);



order_vec_orbit = [0:order_orbit]';
if trig ==1
    Cd_total = repmat(An_total, 1, order_orbit+1).*repmat(cos(order_vec_orbit'*E),numel(An_total),1);
elseif trig ==2
    Cd_total = repmat(An_total, 1, order_orbit+1).*repmat(sin(order_vec_orbit'*E),numel(An_total),1);
elseif trig ==3
    Cd_total = repmat(Bn_total, 1, order_orbit+1).*repmat(cos(order_vec_orbit'*E),numel(An_total),1);
elseif trig ==4
    Cd_total = repmat(Bn_total, 1, order_orbit+1).*repmat(sin(order_vec_orbit'*E),numel(An_total),1);
else
    Cd_total = [An_total,Bn_total];
end
end