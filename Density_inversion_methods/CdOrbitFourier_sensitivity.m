%% Orbit variation of the body-fixed coefficients
clc
clear all
close all
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/DataTools/mice/mice/src/mice/')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/DataTools/mice/mice/lib/' )
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Model 5- 2d fourier series')
%% Initialize constants
Constants_bodyFourier
%% Fourier coefficients
for nn = 1
    Af_total(:,nn) = integral(@(x) cd_environment_exp(x,parameters,1),0, 2*pi,'AbsTol',1e-12,'RelTol',1e-12,'ArrayValued',true);
    Cf_total(:,nn) = integral(@(x) cd_environment_exp(x,parameters,2),0, 2*pi,'AbsTol',1e-12,'RelTol',1e-12,'ArrayValued',true);
end
Af_total_mat = Af_total'/pi;
Cf_total_mat = Cf_total'/pi;
Af_total_mat(1,:) = Af_total_mat(1,:)/2;


Bf_total_mat = zeros(1,Order_o+1);
Df_total_mat = zeros(1,Order_o+1);
Af_std = 0.5*abs(Af_total_mat);
Bf_std = 0.5*abs(Bf_total_mat);
Cf_std = 0.5*abs(Cf_total_mat);
Df_std = 0.5*abs(Df_total_mat);
%%
% angle_E = (0:1:360)*pi/180;
% for nn = 1:numel(Kl_mat)
%     Kl = Kl_mat(nn);
%     for kk = 1:numel(angle_E)
%         [Cd_total(kk),frac(kk),P_o(kk),M(kk),Vi(kk),Ht(kk)] = cd_environment_exp(x,parameters,0);
%         Cd_f(kk) = Af_total_mat(:,nn)'*cos(order_vec*angle_E(kk)) + Bf_total_mat(:,nn)'*sin(order_vec*angle_E(kk));
%         Cd_err(nn,kk) = Cd_total(kk) - Cd_f(kk);
%         
%     end
%     Cdtotal_mat(nn,:) = Cd_total;
%     frac_mat(nn,:) = frac;
%     Po_mat(nn,:) = P_o;
% end
% Cd_err_mean = mean(Cd_err);
% Cd_err_std= std(Cd_err);
%% Plots
% figure(1)
% plot(Cd_err)
% xlabel('Angle')
% ylabel('Cd error')
% set(gca,'FontSize',15)
% title('Drag coefficient error with Fourier series')
