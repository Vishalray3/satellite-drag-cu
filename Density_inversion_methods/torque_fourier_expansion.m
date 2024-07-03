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
Hp_ind = 305;
Ha_ind = 310;
Constants_torque
clear Af_total_mat Bf_total_mat Cf_total_mat Df_total_mat
N_orb = 0;
shape_model = 'plate_dria';           %% 'plate_dria' or 'plate_quasi' or 'sphere'
parameters.shape_model = shape_model;
parameters.N_orb = N_orb;
thresh = 1e-6;
parameters.time_prop = time_prop;
%% Fourier coefficients
frac = [0.1];
Alpha = [1];


parameters.frac = frac;
parameters.Alpha = Alpha;
[Pf_total, Qf_total, An_total, Bn_total] = drag_torque_exp(0,parameters);

Pf_total(abs(Pf_total)< thresh) = 0;
Qf_total(abs(Qf_total)< thresh) = 0;
An_cd = sum(An_total,2);
Bn_cd = sum(Bn_total,2);
An_cd(abs(An_cd)<thresh) = 0;
Bn_cd(abs(Bn_cd)<thresh) = 0;

rho = 2e-12;
ovec = [0:Order_b];
ovec_cd = [0:Order_b+1]';
[r_eci, v_eci] = coe2rv(a_sma,e,inc,raan,w_arg,true_ano,mu_e,str_special, u_arg);
X_state = [r_eci;v_eci];
for n = 1:360
    theta = n-1;

Vi = Fwind*norm(v_eci);
V_vec = [cosd(theta);sind(theta);0];
V_mat = repmat(V_vec, 1, N_plates);
Cd_plates = sum(An_total.*repmat(cosd(ovec_cd*theta),1,N_plates),1) + sum(Bn_total.*repmat(sind(ovec_cd*theta),1,N_plates),1); 
torque_actual(:,n) = 0.5*rho*Vi^2/(2*mass)*sum(cross(Cp, V_mat, 1).*repmat(Cd_plates,3,1) ,2);
torque_fourier(:,n) = 0.5*rho*Vi^2/(2*mass)*(sum(Pf_total.*repmat(cosd(ovec*theta),3,1),2) + ...
    sum(Qf_total.*repmat(sind(ovec*theta),3,1),2));
end

%% 
figure(1)
plot(torque_actual(1,:),'LineWidth',1)
hold on
plot(torque_actual(2,:),'LineWidth',1)
plot(torque_actual(3,:),'LineWidth',1)
plot(torque_fourier(1,:),'--','LineWidth',1)
plot(torque_fourier(2,:),'--','LineWidth',1)
plot(torque_fourier(3,:),'--','LineWidth',1)
legend('True (x)','True (y)','True (z)','Fourier (x)','Fourier (y)','Fourier (z)')
xlabel('Body angle (deg)')
ylabel('Torque (Nm)')
title('Aerodynamic torques')
set(gca,'FontSize',16)
grid on

figure(2)
plot(Pf_total(1,:),'LineWidth',1)
hold on
plot(Pf_total(2,:),'LineWidth',1)
plot(Pf_total(3,:),'LineWidth',1)
plot(Qf_total(1,:),'--','LineWidth',1)
plot(Qf_total(2,:),'--','LineWidth',1)
plot(Qf_total(3,:),'--','LineWidth',1)
legend('Px','Py','Pz','Qx','Qy','Qz')
xlabel('Order')
ylabel('Fourier coefficients')
title('Fourier coefficients')
set(gca,'FontSize',16)
grid on