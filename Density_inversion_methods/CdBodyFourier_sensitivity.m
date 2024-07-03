%% Orbit variation of the body-fixed coefficients
clc
clear all
close all
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/DataTools/mice/mice/src/mice/')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/DataTools/mice/mice/lib/' )
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Model 5- 2d fourier series')
% load nrlmsis_data_2017
%% Initialize constants
Constants_GPM
N = Order_b;
shape_model = 'plate_dria';           %% 'plate_dria' or 'plate_quasi' or 'sphere'   
%% Fourier coefficients
angle_E = (0:0.1:360)*pi/180;

for nn = 1:numel(angle_E)
      
    [Af_analytical(:,nn),Bf_analytical(:,nn),Ht(nn),Vi(nn,1),frac(nn)] = cd_body_exp_analytical(angle_E(nn),e,a_sma,inc,raan,w_arg,mu_e,str_special,...
        u_arg,R,Kl,atm_mass,M_s,Re,Fwind,rho0_all,Hp,H_scale_all,Talt,Ar,Area_plates,Tw,Alpha,N,area_vec,flag_axis);
%                            rho(nn,1) = rho0*exp((Hp - Ht(nn))./H_scale);
                           
end
% Af_avg_rho = Af_analytical*rho/sum(rho);
% Af_avg_full = Af_analytical*(rho.*Vi.^2)/sum(rho.*Vi.^2);
t_lin = (doy-doy(1))*86400 + sod - sod(1);
for kk = 1:numel(time_pred)
    t = t_lin(kk);
    X_state = [reci(:,kk);veci(:,kk)];
    parameters.X_state = X_state;
    parameter.theta = theta(kk);
[rho,rhodiff,Cd_val(kk),BFF_coeff, S, m_r, frac,r_ads] = cd_truth_bff(t,parameters);
end
Cd_ref = Cd_eric.*Aref_eric;
%%
% theta_angle = (0:1:360)*pi/180;
% for nn = 1:numel(angle_E)
%     for kk = 1:numel(theta_angle)
%         [Cd_total(kk),frac(kk),P_o(kk)] = cd_body_exp(angle_E(nn),e,a_sma,inc,raan,w_arg,mu_e,str_special,u_arg,R,Kl...
%             ,shape_model,atm_mass,M_s,Re,Fwind,rho0_all,Hp,H_scale_all,Talt,n_circ,Ar,flag_axis,phi,theta_angle(kk),area_vec,Area_plates,Tw,Alpha,order_vec,0);
%         Cd_f(kk) = Af_total_mat(:,nn)'*cos(order_vec*theta_angle(kk)) + Bf_total_mat(:,nn)'*sin(order_vec*theta_angle(kk));
%         Cd_err(nn,kk) = Cd_total(kk) - Cd_f(kk);
%         
%     end
%     Cdtotal_mat(nn,:) = Cd_total;
%     frac_mat(nn,:) = frac;
%     Po_mat(nn,:) = P_o;
% end
% % Cd_err_mean = mean(Cd_err);
% % Cd_err_std= std(Cd_err);
%% Plots

