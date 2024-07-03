%% Least-squares to calculate K from estimated BFF coefficients
%% Orbit variation of the body-fixed coefficients
clc
clear all
close all
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/DataTools/mice/mice/src/mice/')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/DataTools/mice/mice/lib/' )
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Model 5- 2d fourier series')

% addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Observability study/Results/OFF/Exponential')
%% Initialize constants
%% Fourier coefficients
Constants_bodyFourier

load('Case2_gpm_AllFourierEst','X_nom','P_up','Xf_name_all','mc_Aferr_est','Xf_name','Af_total','Bf_total','Af0_true','Xf_true')

K_true = 1e6;
Order_o = 10;
Order_b = 30;
N = Order_b;                                % highest order of coefficients
order_vec = [0:N]';
order_mat = repmat(order_vec,1,numel(theta));
theta_mat = order_mat.*theta;

N_orb = 10;
bx = a_sma*e/H_scale;
for nn = 1:N_orb+1
    Ino(nn,1) = besseli(nn-1, bx);
end

%% Names
mat_b = repmat([0:Order_b]',1,Order_o+1); mat_o = repmat([0:Order_o],Order_b+1,1);
Af_name = strcat('A', sprintfc('%d',mat_b),sprintfc('%d',mat_o));
Bf_name = strcat('B', sprintfc('%d',mat_b),sprintfc('%d',mat_o));
Cf_name = strcat('C', sprintfc('%d',mat_b),sprintfc('%d',mat_o));
Df_name = strcat('D', sprintfc('%d',mat_b),sprintfc('%d',mat_o));
%% Initial estimates
iter = 10;
thresh = 1e-12;
str_est = [{'A20'},{'B10'}];
ind_est = find(ismember(Xf_name_all,str_est));
R_Af = P_up(10:end,10:end);               %%%%%%%%%%%%% need to do it manually
R_Af = R_Af(ind_est,ind_est);
R_inv = inv(R_Af);

Xf_est_all = [X_nom(10:end)]; % [Af0_true;Xf_true']; %
Xf_est = Xf_est_all(ind_est);
K_old = 1e5;
Kl_est_mat(1) = K_old;

for nn = 1:iter
    Kl = K_old;
    shape_case = 'plate_jac';
    Af_jac = integral(@(x) cd_bodf_exp(x,e,a_sma,inc,raan,w_arg,mu_e,str_special,u_arg,R,Kl,atm_mass,M_s,Re,Fwind,rho0_all,Hp,...
        H_scale_all,Talt,Ar,Area_plates,Tw,Alpha,N,N_orb,area_vec,flag_axis,shape_case,1),0, 2*pi,'AbsTol',1e-12,'RelTol',1e-12,'ArrayValued',true); % cos(body)cos(orbit)
    Cf_jac = integral(@(x) cd_bodf_exp(x,e,a_sma,inc,raan,w_arg,mu_e,str_special,u_arg,R,Kl,atm_mass,M_s,Re,Fwind,rho0_all,Hp,...
        H_scale_all,Talt,Ar,Area_plates,Tw,Alpha,N,N_orb,area_vec,flag_axis,shape_case,2),0, 2*pi,'AbsTol',1e-12,'RelTol',1e-12,'ArrayValued',true); % cos(body)sin(orbit)
    Bf_jac = integral(@(x) cd_bodf_exp(x,e,a_sma,inc,raan,w_arg,mu_e,str_special,u_arg,R,Kl,atm_mass,M_s,Re,Fwind,rho0_all,Hp,...
        H_scale_all,Talt,Ar,Area_plates,Tw,Alpha,N,N_orb,area_vec,flag_axis,shape_case,3),0, 2*pi,'AbsTol',1e-12,'RelTol',1e-12,'ArrayValued',true); % sin(body)cos(orbit)
    Df_jac = integral(@(x) cd_bodf_exp(x,e,a_sma,inc,raan,w_arg,mu_e,str_special,u_arg,R,Kl,atm_mass,M_s,Re,Fwind,rho0_all,Hp,...
        H_scale_all,Talt,Ar,Area_plates,Tw,Alpha,N,N_orb,area_vec,flag_axis,shape_case,4),0, 2*pi,'AbsTol',1e-12,'RelTol',1e-12,'ArrayValued',true); % sin(body)sin(orbit)
    Af_jac = Af_jac/pi;
    Bf_jac = Bf_jac/pi;
    Cf_jac = Cf_jac/pi;
    Df_jac = Df_jac/pi;
    
    Af_jac(:,1) = Af_jac(:,1)/2;
    Bf_jac(:,1) = Bf_jac(:,1)/2;
    Af_ind = ~~Af_jac(1:Order_b+1,1:Order_o+1);
    Bf_ind = ~~Bf_jac(1:Order_b+1,1:Order_o+1);
    Cf_ind = ~~Cf_jac(1:Order_b+1,1:Order_o+1);
    Df_ind = ~~Df_jac(1:Order_b+1,1:Order_o+1);
    Xf_name_all_jac = [Af_name(Af_ind); Bf_name(Bf_ind); Cf_name(Cf_ind); Df_name(Df_ind)];
    Xf_jac = [Af_jac(Af_ind); Bf_jac(Bf_ind); Cf_jac(Cf_ind); Df_jac(Df_ind)];
    %     Af_jac = Af_jac*Ino/Ino(1);
    ind_jac = find(ismember(Xf_name_all_jac,str_est));
    Xf_jac = Xf_jac(ind_jac);
    
    shape_case = 'plate';
    Af_pred = integral(@(x) cd_bodf_exp(x,e,a_sma,inc,raan,w_arg,mu_e,str_special,u_arg,R,Kl,atm_mass,M_s,Re,Fwind,rho0_all,Hp,...
        H_scale_all,Talt,Ar,Area_plates,Tw,Alpha,N,N_orb,area_vec,flag_axis,shape_case,1),0, 2*pi,'AbsTol',1e-12,'RelTol',1e-12,'ArrayValued',true); % cos(body)cos(orbit)
    Cf_pred = integral(@(x) cd_bodf_exp(x,e,a_sma,inc,raan,w_arg,mu_e,str_special,u_arg,R,Kl,atm_mass,M_s,Re,Fwind,rho0_all,Hp,...
        H_scale_all,Talt,Ar,Area_plates,Tw,Alpha,N,N_orb,area_vec,flag_axis,shape_case,2),0, 2*pi,'AbsTol',1e-12,'RelTol',1e-12,'ArrayValued',true); % cos(body)sin(orbit)
    Bf_pred = integral(@(x) cd_bodf_exp(x,e,a_sma,inc,raan,w_arg,mu_e,str_special,u_arg,R,Kl,atm_mass,M_s,Re,Fwind,rho0_all,Hp,...
        H_scale_all,Talt,Ar,Area_plates,Tw,Alpha,N,N_orb,area_vec,flag_axis,shape_case,3),0, 2*pi,'AbsTol',1e-12,'RelTol',1e-12,'ArrayValued',true); % sin(body)cos(orbit)
    Df_pred = integral(@(x) cd_bodf_exp(x,e,a_sma,inc,raan,w_arg,mu_e,str_special,u_arg,R,Kl,atm_mass,M_s,Re,Fwind,rho0_all,Hp,...
        H_scale_all,Talt,Ar,Area_plates,Tw,Alpha,N,N_orb,area_vec,flag_axis,shape_case,4),0, 2*pi,'AbsTol',1e-12,'RelTol',1e-12,'ArrayValued',true); % sin(body)sin(orbit)
    Af_pred = Af_pred/pi;
    Bf_pred = Bf_pred/pi;
    Cf_pred = Cf_pred/pi;
    Df_pred = Df_pred/pi;
    
    Af_pred(:,1) = Af_pred(:,1)/2;
    Bf_pred(:,1) = Bf_pred(:,1)/2;
    Af_ind = ~~Af_pred(1:Order_b+1,1:Order_o+1);
    Bf_ind = ~~Bf_pred(1:Order_b+1,1:Order_o+1);
    Cf_ind = ~~Cf_pred(1:Order_b+1,1:Order_o+1);
    Df_ind = ~~Df_pred(1:Order_b+1,1:Order_o+1);
    Xf_name_all_pred = [Af_name(Af_ind); Bf_name(Bf_ind); Cf_name(Cf_ind); Df_name(Df_ind)];
    Xf_pred = [Af_pred(Af_ind); Bf_pred(Bf_ind); Cf_pred(Cf_ind); Df_pred(Df_ind)];
    %     Af_jac = Af_jac*Ino/Ino(1);
    ind_pred = find(ismember(Xf_name_all_pred,str_est));
    %     Af_pred = Af_pred*Ino/Ino(1);
    Xf_pred = Xf_pred(ind_pred);
    residual = Xf_est - Xf_pred;
    K_var = inv(Xf_jac'*R_inv*Xf_jac);
    K_est = K_old + K_var*Xf_jac'*R_inv*residual;
    if K_est < 0
        K_est = 0;
    end
    K_old = K_est;
    Kl_est_mat(nn+1) = K_old;
end