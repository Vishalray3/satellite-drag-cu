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
Constants_simSpire
clear Af_total_mat Bf_total_mat Cf_total_mat Df_total_mat
N_orb = 0;
shape_model = 'plate_dria';           %% 'plate_dria' or 'plate_quasi' or 'sphere'
parameters.shape_model = shape_model;
parameters.N_orb = N_orb;
thresh = 1e-6;
parameters.time_prop = time_prop;
roll = parameters.roll(1);
pitch = parameters.pitch(1);
yaw = parameters.yaw(1);

order_fourier = 30;
ovec = [0:order_fourier];
%% Fourier coefficients
X_state = X_init;
% ECI TO LVLH
z_hat = -X_state(1:3)/norm(X_state(1:3));   % nadir
y_hat = cross(X_state(4:6), X_state(1:3));  % negative orbit normal
y_hat = y_hat/norm(y_hat);
x_hat = cross(y_hat, z_hat);
x_hat = x_hat/norm(x_hat);

ECI2LVLH = [x_hat,y_hat,z_hat]';

%% LVLH2SBF
roll_mat = [1,0,0; 0,cosd(roll),sind(roll); 0,-sind(roll),cosd(roll)];
pitch_mat = [cosd(pitch), 0, -sind(pitch); 0,1,0; sind(pitch), 0, cosd(pitch)];
yaw_mat = [cosd(yaw), sind(yaw), 0; -sind(yaw),cosd(yaw),0;0,0,1];
LVLH2SBF = roll_mat*pitch_mat*yaw_mat;

ECI2SBF = LVLH2SBF*ECI2LVLH;
SBF2ECI = ECI2SBF';

sun_sbf = ECI2SBF*sun_pos(:,1);
sun_unit = sun_sbf/norm(sun_sbf);
sdel = atan2(sun_unit(2),sqrt(sun_unit(1)^2+sun_unit(3)^2));
slam = atan2(sun_unit(1),sun_unit(3));



Af_total = integral(@(x) srp_nadir_fourier(x,parameters,sdel, order_fourier,1),0, 2*pi,'AbsTol',1e-12,'RelTol',1e-12,'ArrayValued',true);
Bf_total = integral(@(x) srp_nadir_fourier(x,parameters,sdel, order_fourier,2),0, 2*pi,'AbsTol',1e-12,'RelTol',1e-12,'ArrayValued',true);

Af_total_mat(:,:) = Af_total/pi;
Bf_total_mat(:,:) = Bf_total/pi;
Af_total_mat(:,1) = Af_total_mat(:,1)/2;
Bf_total_mat(:,1) = Bf_total_mat(:,1)/2;

Af_total_mat(abs(Af_total_mat)<thresh) = 0;
Bf_total_mat(abs(Bf_total_mat)<thresh) = 0;

%%
for kk = 1
%     angle_step = -141.8;
%     slam = angle_step*pi/180;
    
    X_c = sun_pos(:,1);
    X_rel = X_c - X_state(1:3);                 % third body w.r.t spacecraft
    r_rel = norm(X_rel);
    
    
n_hat = parameters.area_vec;
area = parameters.Area_plates;
rho_spec = parameters.rho_spec;
rho_diff = parameters.rho_diff;

sun_vec = [cos(sdel)*sin(slam);sin(sdel);cos(sdel)*cos(slam)]; %X_rel/r_rel;
sun_vec_mat = repmat(sun_vec,1,numel(area));
cos_theta = dot(sun_vec_mat, n_hat);
%             sp_angle = acosd(cos_theta);
cos_theta(cos_theta<0) = 0;
e_coeff = sum(area.*cos_theta.*(1-rho_spec));
n_coeff_mat = repmat(area.*rho_spec.*cos_theta.^2 + area.*rho_diff.*cos_theta/3,3,1);
n_comp = sum(n_coeff_mat.*n_hat,2);
a_srp(:,kk) = -Psun*AU^2*Cr/(mass*r_rel^2)*(e_coeff*sun_vec + 2*n_comp);

a_srp_fourier(:,kk) = -Psun*AU^2*Cr/(mass*r_rel^2)*(sum(Af_total_mat.*repmat(cos(ovec*slam),3,1),2) + ...
    sum(Bf_total_mat.*repmat(sin(ovec*slam),3,1),2));
end

