function [Af_total_mat,Bf_total_mat] = gsim_estimation_func(Order_b,Order_o,Xf_name_all,X_nom,P_up,parameters,time_prop_curr,oxa,Af_total_mat,Bf_total_mat)
flag_diff = [{'f'},{'s'},{'r_ads'}];
diag_vec = ones(1,numel(flag_diff));
Pk_init = diag(diag_vec);
Pkinv_init = inv(Pk_init);
%% Names
% str_est = [{'A10'},{'A20'},{'A30'},{'A50'},{'A70'},{'B10'}];
% str_est = [{'A10'},{'A20'},{'A30'},{'A40'},{'A50'},{'A60'},{'A70'},{'A80'},{'B10'}];
str_est = [{'A10'},{'A20'},{'A30'},{'A50'},{'A70'},{'B10'}];
mat_b = repmat([0:Order_b]',1,Order_o+1); mat_o = repmat([0:Order_o],Order_b+1,1);
Af_name = strcat('A', sprintfc('%d',mat_b),sprintfc('%d',mat_o));
Bf_name = strcat('B', sprintfc('%d',mat_b),sprintfc('%d',mat_o));
Cf_name = strcat('C', sprintfc('%d',mat_b),sprintfc('%d',mat_o));
Df_name = strcat('D', sprintfc('%d',mat_b),sprintfc('%d',mat_o));

Xf_name_all_jac = [Af_name; Bf_name; Cf_name; Df_name];
ind_name = find(ismember(Xf_name_all_jac,str_est));
%% Initialization
iter = 5;
thresh = 1e-6;
ind_est = find(ismember(Xf_name_all,str_est));
R_Af = P_up(10:end,10:end);               %%%%%%%%%%%%% need to do it manually
R_Af = R_Af(ind_est,ind_est);
R_inv = inv(R_Af);

Xf_est_all = [X_nom(10:end)]; % [Af0_true;Xf_true']; %
Xf_est = Xf_est_all(ind_est);
x0_est = zeros(numel(flag_diff),1);
x_est = zeros(numel(flag_diff),1);
%% Initial Estimates
[~, ~, M, rho_o, T] = density_output(0, parameters.epoch, parameters.doy, time_prop_curr, X_nom(1:6,1), parameters.sun_pos, parameters.eps, parameters);

v_rel = X_nom(4:6,1) - cross([0;0;parameters.omega_e],X_nom(1:3,1));
Vi = norm(v_rel);

Ri = parameters.R/M;
S = Vi./sqrt(2.*Ri.*T);

P_o = rho_o/(oxa*parameters.amu)*parameters.k_b*T;
m_r = M./parameters.M_s;

frac = parameters.Kl.*P_o./(1+parameters.Kl.*P_o);

r_ads = sqrt(0.5*(1+parameters.Alpha*(4*Ri*parameters.Tw./Vi.^2 -1)));
Alpha_s = 3.6*m_r./(1+m_r).^2;
r_s = sqrt(0.5*(1+Alpha_s*(4*Ri*parameters.Tw./Vi.^2 -1)));
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
        [Af_jac, Bf_jac] = bff_analytical(parameters.Area_plates, parameters.Ar, parameters.area_vec, S, r_ads, r_s, frac, Order_b,parameters.flag_axis,...
            flag_method,flag_diff{jj});
        Xf_jac_temp = [Af_jac; Bf_jac];
        Xf_jac_temp = Xf_jac_temp(ind_name);
        Xf_jac = [Xf_jac,Xf_jac_temp];
    end
    
    flag_method = 'coeff';
    [Af_pred, Bf_pred] = bff_analytical(parameters.Area_plates, parameters.Ar, parameters.area_vec, S, r_ads, r_s, frac, Order_b,parameters.flag_axis,...
        flag_method,flag_diff);
    
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
indA = abs(Af_total_mat(:,1))<thresh;
indB = abs(Bf_total_mat(:,1))<thresh;
Af_total_mat = Af_pred;
Bf_total_mat = Bf_pred;
Af_total_mat(indA) = 0;
Bf_total_mat(indB) = 0;