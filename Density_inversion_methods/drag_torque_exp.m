%% function
function [P_total, Q_total, An_total, Bn_total] = drag_torque_exp(angle_step,parameters)
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
Cp = parameters.Cp;

sin_ano = sqrt(1-e^2)*sin(E)./(1-e*cos(E));
cos_ano = (cos(E) - e)./(1-e*cos(E));
true_ano = atan2d(sin_ano,cos_ano);
u_arg = true_ano; % + w_arg
M_ano = E-e*sin(E);
t = M_ano/parameters.n_mean;
time_et = parameters.time_prop(1) + t;  % jd in seconds
[r_eci, v_eci] = coe2rv(a_sma,e,inc,raan,w_arg,true_ano,mu_e,str_special, u_arg);
X_state = [r_eci;v_eci];

% parameters.TEI = time2rotmat_mex(parameters.eop, time_et,X_state, parameters.eqeterms, parameters.nut80);
theta0 = JD2GAST(time_et/86400)*pi/180;
parameters.TEI = [cos(theta0), sin(theta0), 0;...
    -sin(theta0), cos(theta0), 0;...
    0,  0, 1];
[~, ~, M, rho_o, T] = density_output(t, parameters.epoch, parameters.doy, time_et, X_state, parameters.sun_pos, parameters.eps, parameters);

Vi = Fwind*norm(v_eci);

Ri = R/M;
S = Vi./sqrt(2.*Ri.*T);

P_o = rho_o/(oxa*amu)*k_b*T;
m_r = M./M_s;

frac =  parameters.frac; %Kl.*P_o./(1+Kl.*P_o);%
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

An_ads = fourier_plate_dria(A, Ar, S, r_ads,Ci, order_body+1);
An_s = fourier_plate_dria(A, Ar, S, r_s, Ci,order_body+1);

An_plates = frac*An_ads + (1-frac)*An_s;
if strcmp(shape_case, 'plate_jac')
    An_plates = P_o/(1+Kl*P_o)^2*(An_ads-An_s);
end
order_vec_body = [0:order_body+1]';

cos_delta = cos(order_vec_body.*delta);
sin_delta = sin(order_vec_body.*delta);

An_total = An_plates.*cos_delta;
Bn_total = An_plates.*sin_delta;

if strcmp(flag_axis, 'z')
    F_vec = [An_total(2,:); Bn_total(2,:); zeros(1, numel(delta))];
    P_total(:,1)  = sum(cross(Cp, F_vec, 1),2);
    Q_total(:,1) = [0;0;0];
    
    if order_body > 0
        F_vec = [(2*An_total(1,:) + An_total(3,:))/2; Bn_total(3,:)/2; zeros(1, numel(delta))];
        P_total(:,2)  = sum(cross(Cp, F_vec, 1),2);
        F_vec = [Bn_total(3,:)/2; (2*An_total(1,:) - An_total(3,:))/2; zeros(1, numel(delta))];
        Q_total(:,2)  = sum(cross(Cp, F_vec, 1),2);
        if order_body > 1
            for n = 3:order_body+1
                F_vec = [(An_total(n+1,:) + An_total(n-1,:))/2; (Bn_total(n+1,:) - Bn_total(n-1,:))/2; zeros(1, numel(delta))];
                P_total(:,n)  = sum(cross(Cp, F_vec, 1),2);
                F_vec = [(Bn_total(n+1,:) + Bn_total(n-1,:))/2; (An_total(n-1,:) - An_total(n+1,:))/2; zeros(1, numel(delta))];
                Q_total(:,n)  = sum(cross(Cp, F_vec, 1),2);
            end
            
        end
        
    end
    
end

end