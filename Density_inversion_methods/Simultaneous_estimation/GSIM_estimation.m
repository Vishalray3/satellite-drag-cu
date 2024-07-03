%% Least-squares to calculate K from estimated BFF coefficients
%% Orbit variation of the body-fixed coefficients
% clc
% clear all
% close all
% addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/DataTools/mice/mice/src/mice/')
% addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/DataTools/mice/mice/lib/' )
% addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Model 5- 2d fourier series')
% addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Density inversion methods')
% addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Orbit - Processing real data/Full Fourier model/JB08')

% addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Observability study/Results/OFF/Exponential')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Density inversion methods/Simultaneous estimation/spire')
%% Fourier coefficients
clc
clearvars
Constants_spire

% load('1_ekf_msis_cr_cdorder2_iter2','X_nom','Ps_up','Xf_name_all','Xf_name','Af_tseries','Bf_tseries','Af_total_mat','Bf_total_mat')
load('est_msisFourier_srpcball','X_nom','Ps_up','Xf_name_all','Xf_name','Af_total','Bf_total', 'Af_total_mat','Bf_total_mat','X_name_nom','K_ind')

% Order_b = 30;
% Order_o = 0;
N = Order_b;                                % highest order of coefficients
% order_vec = [0:N]';
% order_mat = repmat(order_vec,1,numel(theta));
% theta_mat = order_mat.*theta;


flag_diff = [{'f'}]; %[{'s'},{'r_ads'},{'f'}];
diag_vec = ones(1,numel(flag_diff));
Pk_init = diag(diag_vec);
Pkinv_init = inv(Pk_init);
%% Names
% str_est = [{'A10'},{'A20'},{'A30'},{'A50'},{'A70'},{'B10'}];
% str_est = [{'A10'},{'A20'},{'A30'},{'A40'},{'A50'},{'A60'},{'A70'},{'A80'},{'B10'}];
% str_est = [{'A10'},{'A20'},{'A30'},{'A50'},{'A70'},{'B10'}];
str_est = [{'A10'},{'A20'},{'B10'},{'B20'}];
mat_b = repmat([0:Order_b]',1,Order_o+1); mat_o = repmat([0:Order_o],Order_b+1,1);
Af_name = strcat('A', sprintfc('%d',mat_b),sprintfc('%d',mat_o));
Bf_name = strcat('B', sprintfc('%d',mat_b),sprintfc('%d',mat_o));
Cf_name = strcat('C', sprintfc('%d',mat_b),sprintfc('%d',mat_o));
Df_name = strcat('D', sprintfc('%d',mat_b),sprintfc('%d',mat_o));

Xf_name_all_jac = [Af_name(:); Bf_name(:); Cf_name(:); Df_name(:)];
ind_name = find(ismember(Xf_name_all_jac,str_est));
%% Initialization
iter = 15;
thresh = 1e-6;
ind_est = find(ismember(X_name_nom,str_est));
R_Af = Ps_up(11,11);               %%%%%%%%%%%%% need to do it manually
% R_Af = R_Af(ind_est,ind_est);
R_inv = inv(R_Af);

% if ismember('Cr',estimated_coeff)
% Xf_est_all = [X_nom(11:end)]; % [Af0_true;Xf_true']; %
% else
%     Xf_est_all = [X_nom(10:end)];
% end
Xf_est = X_nom(ind_est);
x0_est = zeros(numel(flag_diff),1);
x_est = zeros(numel(flag_diff),1);
%% Initial Estimates
time_et = time_prop(1);
theta0 = JD2GAST(time_et/86400)*pi/180;
parameters.TEI = [cos(theta0), sin(theta0), 0;...
    -sin(theta0), cos(theta0), 0;...
    0,  0, 1];
[~, ~, M, rho_o, T] = density_output(0, parameters.epoch, parameters.doy, time_prop(1), X_nom(1:6,1), parameters.sun_pos, parameters.eps, parameters);

Vi = norm(X_nom(4:6,1));

Ri = R/M;
S = Vi./sqrt(2.*Ri.*T);

P_o = rho_o/(oxa*amu)*k_b*T;
m_r = M./M_s;

frac = 0.98; %Kl.*P_o./(1+Kl.*P_o);

r_ads = sqrt(0.5*(1+Alpha*(4*Ri*Tw./Vi.^2 -1)));
Alpha_s = 3.6*m_r./(1+m_r).^2;
r_s = sqrt(0.5*(1+Alpha_s*(4*Ri*Tw./Vi.^2 -1)));
% r_s = mean(r_s);

ind_p = find(ismember(flag_diff,'s'));
if ind_p>0
    par_old(ind_p,1) = S;
end
ind_p = find(ismember(flag_diff,'r_ads'));
if ind_p>0
    par_old(ind_p,1) = r_ads;
end
ind_p = find(ismember(flag_diff,'r_s'));
if ind_p>0
    par_old(ind_p,1) = r_s;
end
ind_p = find(ismember(flag_diff,'f'));
if ind_p>0
    par_old(ind_p,1) = frac;
end
par_est_mat(:,1) = par_old;
%% LEAST-SQUARES
for nn = 1:iter
    x0_est = x0_est - x_est;
    ind_p = find(ismember(flag_diff,'s'));
    if ind_p>0
        S = par_old(ind_p);
    end
    ind_p = find(ismember(flag_diff,'r_ads'));
    if ind_p>0
        r_ads = par_old(ind_p);
    end
    ind_p = find(ismember(flag_diff,'r_s'));
    if ind_p>0
        r_s = par_old(ind_p);
    end
    ind_p = find(ismember(flag_diff,'f'));
    if ind_p>0
        frac = par_old(ind_p);
    end
    
    flag_method = 'diff';
    Xf_jac = [];
    for jj = 1:numel(flag_diff)
        [Af_jac, Bf_jac] = bff_analytical(Area_plates, Ar, area_vec, S, r_ads, r_s, frac, Order_b,flag_axis,flag_method,flag_diff{jj});
        Xf_jac_temp = [Af_jac; Bf_jac];
        Xf_jac_temp = Xf_jac_temp(ind_name);
        Xf_jac = [Xf_jac,Xf_jac_temp];
    end
    
    flag_method = 'coeff';
    [Af_pred, Bf_pred] = bff_analytical(Area_plates, Ar, area_vec, S, r_ads, r_s, frac, Order_b,flag_axis,flag_method,flag_diff);
    
    Xf_pred = [Af_pred; Bf_pred];
    Xf_pred = Xf_pred(ind_name);
    residual = Xf_est - Xf_pred;
    par_var = inv(Xf_jac'*R_inv*Xf_jac+ Pkinv_init);
    x_est = par_var*(Xf_jac'*R_inv*residual+Pkinv_init*x0_est);
    par_est = par_old + x_est;
    
    if any(par_est < 0)
        par_est(par_est < 0) = -par_est(par_est < 0);
    end
    par_old = par_est;
    par_est_mat(:,nn+1) = par_old;
end
Af_pred(abs(Af_pred)<thresh) = 0;
Bf_pred(abs(Bf_pred)<thresh) = 0;
% indA = abs(Af_total_mat(:,K_ind))<thresh;
% indB = abs(Bf_total_mat(:,K_ind))<thresh;
% Af_pred(indA) = 0;
% Bf_pred(indB) = 0;
% Af_total_mat(indA) = 0;
% Bf_total_mat(indB) = 0;
% Af_tseries(indA,:) = 0;
% Bf_tseries(indB,:) = 0;
% Af_meanerr = mean(Af_total_mat - Af_tseries, 2);
% Af_rmserr = rms(Af_total_mat - Af_tseries,2);
% Bf_meanerr = mean(Bf_total_mat - Bf_tseries, 2);
% Bf_rmserr = rms(Bf_total_mat - Bf_tseries,2);