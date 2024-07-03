%% function
function [rho,rhodiff,Cd_val,Aref,Cl_val,Cd_ads, Cd_s] = cd_truth_bff(time_jd, parameters)
X_state = parameters.X_state;
R = parameters.R;
Kl = parameters.Kl;
M_s = parameters.M_s;
Ar = parameters.Ar;
A = parameters.Area_plates;
Tw = parameters.Tw;
Alpha = parameters.Alpha;
order_body = parameters.Order_b;
Ai = parameters.area_vec;
flag_axis = parameters.flag_axis;
shape_case = parameters.shape_model;
amu = parameters.amu;
oxa = parameters.atm_mass(2);
k_b = parameters.k_b;
phi = parameters.phi;
theta = parameters.theta;
omega_e = parameters.omega_e;
r_eci = X_state(1:3);
v_eci = X_state(4:6);
v_sbf = parameters.v_sbf;

[rho, rhodiff, M, rho_o, T] = density_output(time_jd, X_state, parameters.sun_pos, parameters.eps, parameters);

v_rel = v_eci - cross([0;0;omega_e],r_eci);
Vi = norm(v_rel);

Ri = R/M;
S = Vi./sqrt(2.*Ri.*T);

P_o = rho_o/(oxa*amu)*k_b*T;
m_r = M./M_s;

frac = parameters.frac; %Kl.*P_o./(1+Kl.*P_o); 
if strcmp(shape_case, 'plate_quasi')
    [Cd_val,Cl_val,Aref,Cd_ads,Cd_s] = single_plate_quasi_cd(phi, theta, S, Ai, m_r, Ri, Tw, Vi, A, Ar, Alpha, frac, flag_axis,v_sbf);
elseif strcmp(shape_case, 'plate_dria')
    [Cd_val,Cl_val,Aref,Cd_ads,Cd_s] = single_plate_dria_cd(phi, theta, S, Ai, m_r, Ri, Tw, Vi, A, Ar, Alpha, frac, flag_axis, v_sbf);
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
   elseif strcmp(flag_axis, 'x')
        Ci = sqrt(Ai(2,:).^2+Ai(3,:).^2);
        delta = atan2(Ai(2,:), Ai(3,:));          
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
    
    An_total = sum(An_plates.*cos_delta,2);
    Bn_total = sum(An_plates.*sin_delta,2);
    
    
    BFF_coeff = [An_total,Bn_total];
elseif strcmp(shape_case, 'plate_dria_multi')
    [Cd_val,Aref] = multi_plate_dria_cd(phi, theta, S, Ai, m_r, Ri, Tw, Vi, A, Ar, Alpha, frac, flag_axis, v_sbf);
    Cl_val = 0;
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
    elseif strcmp(flag_axis, 'x')
        Ci = sqrt(Ai(2,:).^2+Ai(3,:).^2);
        delta = atan2(Ai(2,:), Ai(3,:));    
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
    
    An_total = sum(An_plates.*cos_delta,2);
    Bn_total = sum(An_plates.*sin_delta,2);
    
    
    BFF_coeff = [An_total,Bn_total];    
elseif strcmp(shape_case,'sphere') || strcmp(shape_case,'sphere_jac')
    Cd_val = sphere_dria_cd(S, m_r(1), Ri, Tw, Vi, Alpha,frac,Kl,P_o,shape_case);
    BFF_coeff = [0,0];
end



end