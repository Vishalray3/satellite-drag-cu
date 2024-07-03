%% Least-squares to calculate K from estimated OFF coefficients
clc
clear all
close all
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/DataTools/mice/mice/src/mice/')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/DataTools/mice/mice/lib/' )
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Model 5- 2d fourier series')

addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Observability study/Results/OFF/Exponential')
%% Initialize constants
Constants_Fourier

load('Case1_sphere_order3','X_nom','P_up','Af_true','K_true')
N = numel(X_nom(7:end))-1;                                % highest order of coefficients
Af_est = X_nom(7:7+N); % Af_true'; %
R_Af = P_up(7:7+N,7:7+N);
% R_Af = eye(N+1,N+1);
R_inv = inv(R_Af);
order_vec = [0:N]';
order_mat = repmat(order_vec,1,numel(theta));
theta_mat = order_mat.*theta;
%% Initial estimates
K_old = 1e7;
iter = 10;
%% Fourier coefficients
Kl_mat(1) = K_old;
for nn = 1:iter
    nn
    Kl = K_old;
    Af_jac = integral(@(x) cd_environment_exp(x,e,a_sma,inc,raan,w_arg,mu_e,str_special,u_arg,R,Kl,'sphere_jac',atm_mass, M_s,Re,Fwind,rho0_all,Hp...
        ,H_scale_all,Talt,n_circ,Ar,flag_axis,phi,theta,area_vec,Area_plates,Tw,Alpha,order_vec,1),0, 2*pi,'AbsTol',1e-12,'RelTol',1e-12,'ArrayValued',true);
    Af_jac = Af_jac/pi;
    Af_jac(1) = Af_jac(1)/2;

    Af_pred = integral(@(x) cd_environment_exp(x,e,a_sma,inc,raan,w_arg,mu_e,str_special,u_arg,R,Kl,'sphere',atm_mass, M_s,Re,Fwind,rho0_all,Hp...
        ,H_scale_all,Talt,n_circ,Ar,flag_axis,phi,theta,area_vec,Area_plates,Tw,Alpha,order_vec,1),0, 2*pi,'AbsTol',1e-12,'RelTol',1e-12,'ArrayValued',true);
    Af_pred = Af_pred/pi;
    Af_pred(1) = Af_pred(1)/2;    
    residual = Af_est - Af_pred
    K_var = inv(Af_jac'*R_inv*Af_jac);
    K_est = K_old + K_var*Af_jac'*R_inv*residual;
    if K_est < 0
        K_est = 0;
    end
    K_old = K_est;
    Kl_mat(nn+1) = K_old;
end
K_error = K_true - K_est;
%%
angle_E = (0:1:360)*pi/180;
for nn = 1:numel(Kl_mat)
    Kl = Kl_mat(nn);
    for kk = 1:numel(angle_E)
        [Cd_total(kk),frac(kk),P_o(kk)] = cd_environment_exp(angle_E(kk),e,a_sma,inc,raan,w_arg,mu_e,str_special,u_arg,R,Kl...
            ,'sphere',atm_mass,M_s,Re,Fwind,rho0_all,Hp,H_scale_all,Talt,n_circ,Ar,flag_axis,phi,theta,area_vec,Area_plates,Tw,Alpha,order_vec,0);
        Cd_f(kk) = Af_est'*cos(order_vec*angle_E(kk));
        Cd_err(nn,kk) = Cd_total(kk) - Cd_f(kk);
        
    end
    Cdtotal_mat(nn,:) = Cd_total;
    frac_mat(nn,:) = frac;
    Po_mat(nn,:) = P_o;
end

%% Plots
figure(1)
plot(Cd_err')
xlabel('Angle')
ylabel('Cd error')
set(gca,'FontSize',15)
title('Drag coefficient error with Fourier series')
