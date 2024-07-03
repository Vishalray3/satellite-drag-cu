%% Orbit variation of the body-fixed coefficients
clc
clear all
close all
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/DataTools/mice/mice/src/mice/')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/DataTools/mice/mice/lib/' )
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Model 5- 2d fourier series')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Orbit - Processing real data/Full Fourier model/JB08')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Density inversion methods/Simultaneous estimation')
addpath ('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/public/quaternions_example_Spire_shaylah')
%% Initialize constants
% Constants_bodyFourier
Constants_satShape
clear Af_total_mat Bf_total_mat Cf_total_mat Df_total_mat
N_orb = 0;
shape_model = 'plate_dria';           %% 'plate_dria' or 'plate_quasi' or 'sphere'
parameters.shape_model = shape_model;
parameters.N_orb = N_orb;
thresh = 1e-6;
parameters.time_prop = time_prop;
%% Fourier coefficients
frac = [1];
Alpha = [1];

% Instantaneous
for nn = 1:numel(frac)
parameters.frac = frac(nn);
parameters.Alpha = Alpha(nn);
Fourier_coeff = cd_bodf_exp(0,parameters,0);
Af_total_mat_inst(:,nn) = Fourier_coeff(1:Order_b+1);
Bf_total_mat_inst(:,nn) = Fourier_coeff(Order_b+2:end);
end

% Orbit averaged
for nn = 1:numel(frac)
    nn
    parameters.frac = frac(nn);
    parameters.Alpha = Alpha(nn);
    Af_total = integral(@(x) cd_bodf_exp(x,parameters,1),0, 2*pi,'AbsTol',1e-12,'RelTol',1e-12,'ArrayValued',true); % cos(body)cos(orbit)
    Cf_total = integral(@(x) cd_bodf_exp(x,parameters,2),0, 2*pi,'AbsTol',1e-12,'RelTol',1e-12,'ArrayValued',true); % cos(body)sin(orbit)
    Bf_total = integral(@(x) cd_bodf_exp(x,parameters,3),0, 2*pi,'AbsTol',1e-12,'RelTol',1e-12,'ArrayValued',true); % sin(body)cos(orbit)
    Df_total = integral(@(x) cd_bodf_exp(x,parameters,4),0, 2*pi,'AbsTol',1e-12,'RelTol',1e-12,'ArrayValued',true); % sin(body)sin(orbit)
    
    Af_total_mat(:,:,nn) = Af_total/pi;
    Cf_total_mat(:,:,nn) = Cf_total/pi;
    Bf_total_mat(:,:,nn) = Bf_total/pi;
    Df_total_mat(:,:,nn) = Df_total/pi;
    Af_total_mat(:,1,nn) = Af_total_mat(:,1,nn)/2;
    Bf_total_mat(:,1,nn) = Bf_total_mat(:,1,nn)/2;
%     
%     Af_Kl(:,nn) = Af_total_mat(:,:,nn)*Ino/Ino(1);
%     Bf_Kl(:,nn) = Bf_total_mat(:,:,nn)*Ino/Ino(1);
end
Af_total_mat(abs(Af_total_mat)<thresh) = 0;
Bf_total_mat(abs(Bf_total_mat)<thresh) = 0;
Cf_total_mat(abs(Cf_total_mat)<thresh) = 0;
Df_total_mat(abs(Df_total_mat)<thresh) = 0;




% Af_Kl(abs(Af_Kl)<thresh) = 0;
% Bf_Kl(abs(Bf_Kl)<thresh) = 0;
% Af_mean = mean(Af_Kl,2)';
% Af_std = std(Af_total_mat,0,'all');%10.^(log10(abs(Af_mean)));
% Bf_mean = mean(Bf_Kl,2)';
% Bf_std = 10.^(log10(abs(Bf_mean)));
% Af_E0 = cd_bodf_exp(0,e,a_sma,inc,raan,w_arg,mu_e,str_special,u_arg,R,Kl,atm_mass,M_s,Re,Fwind,rho0_all,Hp,...
%     H_scale_all,Talt,Ar,Area_plates,Tw,Alpha,N,N_orb,area_vec,flag_axis,0);
% Af_std = 0.5*abs(Af_total_mat);
% Bf_std = 0.5*abs(Bf_total_mat);
% Cf_std = 0.5*abs(Cf_total_mat);
% Df_std = 0.5*abs(Df_total_mat);

% Af_std = abs(Af_total_mat(:,:,5) - Af_total_mat(:,:,4)); 
% Bf_std = abs(Bf_total_mat(:,:,5) - Bf_total_mat(:,:,4)); 
% Cf_std = abs(Cf_total_mat(:,:,5) - Cf_total_mat(:,:,4)); 
% Df_std = abs(Df_total_mat(:,:,5) - Df_total_mat(:,:,4)); 
% Af_std(Af_std<1e-5) = 0;
% Bf_std(Bf_std<1e-5) = 0;
% Cf_std(Cf_std<1e-5) = 0;
% Df_std(Df_std<1e-5) = 0;
%%
% load('spire83_propagation_2018_11_7','ode_step', 'X1','sod_uni')
% % load('spire83_11_7_bff_avg')
% X_true_aug = interp1(ode_step, X1,sod_uni,'spline');
% X_true_aug = X_true_aug';
% % theta_new = smooth(theta, 500,'moving')';
% % phi_new = smooth(phi, 500,'moving')';
% t_lin = time_prop_utc; %(doy-doy(1))*86400 + sod - sod(1);
% theta = repmat(theta, N_plates,1);
% phi = repmat(phi, N_plates,1);
% for nn = 1:numel(frac)
%     parameters.frac = frac(nn);
%     parameters.Alpha = Alpha(nn);
% for kk = 1:30000; %numel(time_prop)
%     rot_ECI2ECEF = time2rotmat_mex(eop, time_prop(kk),X_true_aug(1:6,kk), eqeterms, nut80);
% % rot_ECI2ECEF = cspice_pxform( 'J2000', 'ITRF93', time_var );
%     parameters.TEI = rot_ECI2ECEF;
%     t = t_lin(kk);
%     X_state = X_true_aug(1:6,kk);
%     parameters.X_state = X_state;
%     parameters.theta = pi/180*theta(:,kk)';
%     parameters.phi = pi/180*phi(:,kk)';
%     parameters.time_jd = time_prop(kk);
% [rho,rhodiff,Cd_val(nn,kk),BFF_coeff, S, m_r, ~,r_ads] = cd_truth_bff(t,parameters);
% Cd_f(nn,kk) = Af_total_mat(:,:,nn)'*cosd([0:Order_b]'*theta(6,kk)) + Bf_total_mat(:,:,nn)'*sind([0:Order_b]'*theta(6,kk));
% end
% end
% % Cd_ref = Cd_eric.*Aref_eric;
% 
%% PLots
% load spire83_11_7_bff_avg_msis
% Cd_val1 = Cd_val;
% err_msis = (Cd_val - Cd_f)./Cd_val*100;

% figure(1)
% plot(sod_uni(1:30000)/60, err_msis(1,:), 'LineWidth',1)
% hold on
% plot(sod_uni(1:30000)/60, err_msis(2,:), 'LineWidth',1)
% plot(sod_uni(1:30000)/60, err_msis(3,:), 'LineWidth',1)
% % plot(sod_uni(1:30000)/60, err_msis(4,:), 'LineWidth',1)
% xlabel('Time (mins)')
% ylabel('Relative error (%)')
% leg = legend('DRIA: $f=1, \alpha = 1$', 'DRIA: $f=1, \alpha = 0$', 'DRIA: $f=0, \alpha = \alpha_s$');
% set(leg, 'interpreter','latex','FontSize',15)
% set(gca,'FontSize',14)
% grid on
% title('Drag-coefficient error due to averaging')
% set(gca,'FontSize',14)
% 
% load spire83_11_7_bff_avg_jb08
% err_msis = (Cd_val1 - Cd_val)./Cd_val1*100;
% figure(2)
% plot(sod_uni(1:30000)/60, err_msis(1,:), 'LineWidth',1)
% hold on
% plot(sod_uni(1:30000)/60, err_msis(2,:), 'LineWidth',1)
% plot(sod_uni(1:30000)/60, err_msis(3,:), 'LineWidth',1)
% % plot(sod_uni(1:30000)/60, err_msis(4,:), 'LineWidth',1)
% xlabel('Time (mins)')
% ylabel('Relative error (%)')
% leg = legend('DRIA: $f=1, \alpha = 1$', 'DRIA: $f=1, \alpha = 0$', 'DRIA: $f=0, \alpha = \alpha_s$');
% set(leg, 'interpreter','latex','FontSize',14)
% set(gca,'FontSize',14)
% grid on
% title('Drag-coefficient error between NRLMSISE-00 and JB08')
% set(gca,'FontSize',14)
