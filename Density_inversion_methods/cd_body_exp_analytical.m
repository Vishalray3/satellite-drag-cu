%% function
function [An_total,Bn_total,Ht,Vi,frac] = cd_body_exp_analytical(angle_step,e,a_sma,inc,raan,w_arg,mu_e,str_special,u_arg,R,Kl...
                               ,atm_mass,M_s,Re,Fwind,rho0_all,H0,H_scale_all,T,Ar,A,Tw,Alpha,order,Ai,flag_axis)
global amu oxa k_b
E = angle_step;
sin_ano = sqrt(1-e^2)*sin(E)./(1-e*cos(E));
cos_ano = (cos(E) - e)./(1-e*cos(E));
true_ano = atan2d(sin_ano,cos_ano);
M_ano = E-e*sin(E);
[r_eci, v_eci] = coe2rv(a_sma,e,inc,raan,w_arg,true_ano,mu_e,str_special, u_arg);

Ht = norm(r_eci) - Re;

Vi = Fwind*norm(v_eci);

n_den = rho0_all.*exp((H0 - Ht)./H_scale_all);
M = 1/(n_den/sum(n_den)*(1./atm_mass'));
Ri = R/M;
S = Vi./sqrt(2.*Ri.*T);

rho_o = n_den(:,2) + n_den(:,8);
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
end

An_ads = fourier_plate_dria(A, Ar, S, r_ads,Ci, order);
An_s = fourier_plate_dria(A, Ar, S, r_s, Ci,order);

An_plates = frac*An_ads + (1-frac)*An_s;
order_vec = [0:order]';

cos_delta = cos(order_vec.*delta);
sin_delta = sin(order_vec.*delta); 

An_total = sum(An_plates.*cos_delta,2);
Bn_total = sum(An_plates.*sin_delta,2);
end