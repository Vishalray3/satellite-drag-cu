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
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Density inversion methods/Simultaneous estimation/Torque estimation')
%% Fourier coefficients
clc
clearvars
Hp_ind = 305;
Ha_ind = 310;
Constants_torque

% load('1_ekf_msis_cr_cdorder2_iter2','X_nom','Ps_up','Xf_name_all','Xf_name','Af_tseries','Bf_tseries','Af_total_mat','Bf_total_mat')
load('300km_torquesimulestJB08_iter1','Ps_up','Xf_name1','Xf_name2','Xf_name3','X_est1','X_est2','X_est3','K_ind', 'reci','veci')
X_nom = [reci(:,1);veci(:,1)];

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
str_est1 = [{'Q11'},{'Q13'}];
str_est2 = [{'P21'},{'P23'}];
str_est3 = [{'P31'},{'P32'},{'P33'}];


mat_b = repmat([1:3]',1,Order_b+1); mat_o = repmat([0:Order_b],3,1);
Pf_name = strcat('P', sprintfc('%d',mat_b),sprintfc('%d',mat_o));
Qf_name = strcat('Q', sprintfc('%d',mat_b),sprintfc('%d',mat_o));
Xf1_name = [Pf_name(1,:)'; Qf_name(1,:)']; Xf2_name = [Pf_name(2,:)'; Qf_name(2,:)']; Xf3_name = [Pf_name(3,:)'; Qf_name(3,:)';];
% Xf_name1 = Xf1_name(Xf_ind1);    Xf_name2 = Xf2_name(Xf_ind2);    Xf_name3 = Xf3_name(Xf_ind3);

vec_est1 = find(ismember(Xf1_name,str_est1));    N_est1 = numel(vec_est1);
vec_est2 = find(ismember(Xf2_name,str_est2));    N_est2 = numel(vec_est2);
vec_est3 = find(ismember(Xf3_name,str_est3));    N_est3 = numel(vec_est3);

%% Initialization
iter = 15;
thresh = 1e-6;
ind_est1 = find(ismember(Xf_name1,str_est1));
ind_est2 = find(ismember(Xf_name2,str_est2));
ind_est3 = find(ismember(Xf_name3,str_est3));
P_f = Ps_up(4:10,4:10);               %%%%%%%%%%%%% need to do it manually
R_Af = P_f;
% R_Af = R_Af(ind_est,ind_est);
R_inv = inv(R_Af);

% if ismember('Cr',estimated_coeff)
% Xf_est_all = [X_nom(11:end)]; % [Af0_true;Xf_true']; %
% else
%     Xf_est_all = [X_nom(10:end)];
% end
Xf_est1 = X_est1(ind_est1);
Xf_est2 = X_est2(ind_est2);
Xf_est3 = X_est3(ind_est3);
Xf_est = [Xf_est1; Xf_est2; Xf_est3];
x0_est = zeros(numel(flag_diff),1);
x_est = zeros(numel(flag_diff),1);
%% Initial Estimates
time_et = time_prop(1);
theta0 = JD2GAST(time_et/86400)*pi/180;
parameters.TEI = [cos(theta0), sin(theta0), 0;...
    -sin(theta0), cos(theta0), 0;...
    0,  0, 1];
[~, ~, M, rho_o, T] = density_output(0, parameters.epoch, parameters.doy, time_prop(1), X_nom, parameters.sun_pos, parameters.eps, parameters);

Vi = norm(X_nom(4:6,1));

Ri = R/M;
parameters.S = Vi./sqrt(2.*Ri.*T);


m_r = M./M_s;

parameters.frac = 0.4; %Kl.*P_o./(1+Kl.*P_o);

parameters.r_ads = sqrt(0.5*(1+Alpha*(4*Ri*Tw./Vi.^2 -1)));
Alpha_s = 3.6*m_r./(1+m_r).^2;
parameters.r_s = sqrt(0.5*(1+Alpha_s*(4*Ri*Tw./Vi.^2 -1)));
% r_s = mean(r_s);

ind_p = find(ismember(flag_diff,'s'));
if ind_p>0
    par_old(ind_p,1) = parameters.S;
end
ind_p = find(ismember(flag_diff,'r_ads'));
if ind_p>0
    par_old(ind_p,1) = parameters.r_ads;
end
ind_p = find(ismember(flag_diff,'r_s'));
if ind_p>0
    par_old(ind_p,1) = parameters.r_s;
end
ind_p = find(ismember(flag_diff,'f'));
if ind_p>0
    par_old(ind_p,1) = parameters.frac;
end
par_est_mat(:,1) = par_old;
%% LEAST-SQUARES
for nn = 1:iter
    x0_est = x0_est - x_est;
    ind_p = find(ismember(flag_diff,'s'));
    if ind_p>0
        parameters.S = par_old(ind_p);
    end
    ind_p = find(ismember(flag_diff,'r_ads'));
    if ind_p>0
        parameters.r_ads = par_old(ind_p);
    end
    ind_p = find(ismember(flag_diff,'r_s'));
    if ind_p>0
        parameters.r_s = par_old(ind_p);
    end
    ind_p = find(ismember(flag_diff,'f'));
    if ind_p>0
        parameters.frac = par_old(ind_p);
    end
    
    flag_method = 'diff';
    Xf_jac = [];
    for jj = 1:numel(flag_diff)
        [Pf_jac, Qf_jac] = torque_gsim_inversion(parameters,flag_method,flag_diff{jj});
        Xf_temp1 = [Pf_jac(1,:)'; Qf_jac(1,:)'];
        Xf_temp2 = [Pf_jac(2,:)'; Qf_jac(2,:)'];
        Xf_temp3 = [Pf_jac(3,:)'; Qf_jac(3,:)'];
        
        Xf_jac_temp = [Xf_temp1(vec_est1); Xf_temp2(vec_est2); Xf_temp3(vec_est3)];
        Xf_jac = [Xf_jac,Xf_jac_temp];
    end
    
    flag_method = 'coeff';
    [Pf_pred, Qf_pred] = torque_gsim_inversion(parameters,flag_method,flag_diff);
    
    Xf_pred1 = [Pf_pred(1,:)'; Qf_pred(1,:)'];
    Xf_pred2 = [Pf_pred(2,:)'; Qf_pred(2,:)'];
    Xf_pred3 = [Pf_pred(3,:)'; Qf_pred(3,:)'];
    
    
    Xf_pred = [Xf_pred1(vec_est1); Xf_pred2(vec_est2); Xf_pred3(vec_est3)];
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
Pf_pred(abs(Pf_pred)<thresh) = 0;
Qf_pred(abs(Qf_pred)<thresh) = 0;
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