%% Orbit variation of the body-fixed coefficients
rng('default')
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
N_cases = 500;

%% Errors (relative standard deviations for Gaussan noise) 
del_V = 0.1/3;
del_alpha = 0.8/3;
del_M = 0.5/3;
del_T = 0.5/3;
%% Fourier coefficients
frac = [1];
Alpha = [1];
parameters.frac = frac;
parameters.Alpha = Alpha;
% Instantaneous
for nn = 1:N_cases
random_array = normrnd([0,0,0,0],[del_V, del_alpha, del_M, del_T]);
% random_array = rand(1,4).*[del_V, del_alpha, del_M, del_T]*3;
parameters.errors = random_array;
Fourier_coeff = cd_bodf_exp_mc(0,parameters,0);
Af_total_mat(:,nn) = Fourier_coeff(1:Order_b+1);
Bf_total_mat(:,nn) = Fourier_coeff(Order_b+2:end);
end
Af_stat(:,1) = mean(Af_total_mat,2);
Af_stat(:,2) = std(Af_total_mat,1,2);
Af_stat(:,3) = max(abs(Af_total_mat - Af_stat(:,1)),[],2);
Bf_stat(:,1) = mean(Bf_total_mat,2);
Bf_stat(:,2) = std(Bf_total_mat,1,2);
Bf_stat(:,3) = max(abs(Bf_total_mat - Bf_stat(:,1)),[],2);


% 
% Af_total_mat(abs(Af_total_mat)<thresh) = 0;
% Bf_total_mat(abs(Bf_total_mat)<thresh) = 0;
% Cf_total_mat(abs(Cf_total_mat)<thresh) = 0;
% Df_total_mat(abs(Df_total_mat)<thresh) = 0;



