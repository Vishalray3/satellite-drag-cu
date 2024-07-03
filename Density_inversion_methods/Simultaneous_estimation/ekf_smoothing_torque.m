%% Torque Fourier coefficient estimation
rng('default')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Density inversion methods')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Density inversion methods/Simultaneous estimation/Torque estimation')
clc
clearvars
Hp_ind = 305;
Ha_ind = 310;
Constants_torque
load('torque_fourier_case1_frac_JB08','Pf_total','Qf_total')
Pf_total_err = Pf_total;
Qf_total_err = Qf_total;
% load('gsim_torque_est_iter1','Pf_pred','Qf_pred')
% Pf_total_err = Pf_pred;
% Qf_total_err = Qf_pred;

load('torque_fourier_case1','Pf_total','Qf_total','An_total','Bn_total')
load('300km_torque_case1', 'reci', 'veci','rho','theta','time_prop_utc','omega_e','rho_msis')
load('300km_torquefourier_noise')
theta = theta*180/pi;
Order_b = 30;
Pf_std = 10*abs(Pf_total - Pf_total_err); %0.1*abs(Pf_total);
Qf_std = 10*abs(Qf_total - Qf_total_err); %0.1*abs(Qf_total);
R_aug = diag([1e-14,1e-14,1e-14]);
rho_nom = rho_msis;
% Pf_total_err = Pf_total+ Pf_noise;
% Qf_total_err = Qf_total+ Qf_noise;
%% nominal entries corresponing to the indices
X_f1 = [Pf_total_err(1,:)'; Qf_total_err(1,:)']; X_f2 = [Pf_total_err(2,:)'; Qf_total_err(2,:)']; X_f3 = [Pf_total_err(3,:)'; Qf_total_err(3,:)';];
Xf_ind1 = ~~X_f1; Xf_ind2 = ~~X_f2; Xf_ind3 = ~~X_f3;
X_f1 = X_f1(Xf_ind1); X_f2 = X_f2(Xf_ind2); X_f3 = X_f3(Xf_ind3);
Xf_est = [X_f1;X_f2;X_f3];
N_f1 = numel(X_f1); N_f2 = numel(X_f2); N_f3 = numel(X_f3);

X_f1t = [Pf_total(1,:)'; Qf_total(1,:)']; X_f2t = [Pf_total(2,:)'; Qf_total(2,:)']; X_f3t = [Pf_total(3,:)'; Qf_total(3,:)';];
X_f1t = X_f1t(Xf_ind1); X_f2t = X_f2t(Xf_ind2); X_f3t = X_f3t(Xf_ind3);
Xf_true = [X_f1t;X_f2t;X_f3t];
%% nominal standard deviations
Xf_std1 = [Pf_std(1,:)'; Qf_std(1,:)']; Xf_std2 = [Pf_std(2,:)'; Qf_std(2,:)']; Xf_std3 = [Pf_std(3,:)'; Qf_std(3,:)';];
Xf_std1 = Xf_std1(Xf_ind1);            Xf_std2 = Xf_std2(Xf_ind2);          Xf_std3 = Xf_std3(Xf_ind3);

% P_init_inv = inv(P_init);
%% Names of the coefficients
mat_b = repmat([1:3]',1,Order_b+1); mat_o = repmat([0:Order_b],3,1);
Pf_name = strcat('P', sprintfc('%d',mat_b),sprintfc('%d',mat_o));
Qf_name = strcat('Q', sprintfc('%d',mat_b),sprintfc('%d',mat_o));
Xf1_name = [Pf_name(1,:)'; Qf_name(1,:)']; Xf2_name = [Pf_name(2,:)'; Qf_name(2,:)']; Xf3_name = [Pf_name(3,:)'; Qf_name(3,:)';];
Xf_name1 = Xf1_name(Xf_ind1);    Xf_name2 = Xf2_name(Xf_ind2);    Xf_name3 = Xf3_name(Xf_ind3);
% vectors from those entries, : operator ensures column vector

%% Estimated vector
% str_est1 = [{'P01'},{'Q11'},{'Q13'},{'Q15'}];
% str_est2 = [{'P21'},{'P23'},{'P25'}];
% str_est3 = [{'P31'},{'P32'},{'P33'},{'P35'},{'P37'},{'Q31'}];

% str_est1 = [{'Q11'},{'Q13'}];
% str_est2 = [{'P21'},{'P23'}];
% str_est3 = [{'P31'},{'P32'},{'P33'},{'Q32'}];
str_est1 = [{}];
str_est2 = [{}];
str_est3 = [{}];

% str_est1 = [{'P12'},{'Q11'},{'Q13'}];
% str_est2 = [{'P21'},{'P23'}];
% str_est3 = [{'P31'},{'P32'},{'P33'},{'Q31'},{'Q32'}];
vec_est1 = find(ismember(Xf_name1,str_est1));    N_est1 = numel(vec_est1);
vec_est2 = find(ismember(Xf_name2,str_est2));    N_est2 = numel(vec_est2);
vec_est3 = find(ismember(Xf_name3,str_est3));    N_est3 = numel(vec_est3);

N_f = N_est1+N_est2+N_est3;
P_f = diag([Xf_std1(vec_est1).^2;Xf_std2(vec_est2).^2;Xf_std3(vec_est3).^2]);
Hmat = [];
Ymat = [];
x_prior = zeros(N_f,1);
x_est = x_prior;
N_iter = 7;
N_plates = numel(An_total(1,:));
ovec_cd = [0:Order_b+1]';
ovec = [0:Order_b];

err_init = [X_f1t(vec_est1); X_f2t(vec_est2); X_f3t(vec_est3)] - [X_f1(vec_est1); X_f2(vec_est2); X_f3(vec_est3)];

N_st = 3;
N_est = 3;
vec_est = 1:N_est;
N_all = N_f+N_est;
x_prior_fixed = zeros(N_all,1);
vec_con = [];

P_rho = diag(Xs_std.^2);
P_prior = blkdiag(P_rho(vec_est,vec_est),P_f);
X_nom = X_s;
X_nom_f = [X_f1(vec_est1); X_f2(vec_est2); X_f3(vec_est3)];

stm_init = eye(N_est);
stm_q_init = zeros(N_est);
stm_vec_init = stm_init(:);
stm_q_vec_init = stm_q_init(:);
estimated_coeff.rho_DMC =1;
parameters.N_st = N_st;
parameters.vec_est = vec_est;
%% ODE options
del_T = 10;
ode4_delT = 10;
Tsamp = del_T/ode4_delT;
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
step = time_prop_utc;
ode_step = step(1):ode4_delT:step(end);
%% ekf 

for iter = 1:N_iter
    iteration = iter-1
    
    P_up = P_prior;
    sigma_x(:,1) = sqrt(diag(P_up));
    
    for ii = 1:numel(time_prop_utc)
        V_vec = [cosd(theta(ii)); sind(theta(ii)); 0];
        V_mat = repmat(V_vec, 1, N_plates);
        Cd_plates = sum(An_total.*repmat(cosd(ovec_cd*theta(ii)),1,N_plates),1) + sum(Bn_total.*repmat(sind(ovec_cd*theta(ii)),1,N_plates),1);
        Cd_total(ii) = sum(Cd_plates);
        v_r = veci(:,ii) - cross([0;0;omega_e],reci(:,ii));
        Vi = norm(v_r);
        
        y_meas(:,ii) = -0.5*rho(ii)*Vi^2/(mass)*sum(cross(Cp, V_mat, 1).*repmat(Cd_plates,3,1) ,2);
        
        Xf_arg = [cosd(ovec*theta(ii))';sind(ovec*theta(ii))'];
        Xf_arg1 = Xf_arg(Xf_ind1);   Xf_arg2 = Xf_arg(Xf_ind2);    Xf_arg3 = Xf_arg(Xf_ind3);
        
        X_est1 = Xf_est(1:N_f1); X_est2 = Xf_est(N_f1+1:N_f1+N_f2); X_est3 = Xf_est(N_f1+N_f2+1:end);
        Tx = sum(X_est1.*Xf_arg1); Ty = sum(X_est2.*Xf_arg2); Tz = sum(X_est3.*Xf_arg3);
        
        
        
        if time_prop_utc_ekf(ii) == 0
            k = 0;
            X_aug(:,1) = [X_nom(:,1);  stm_vec_init; stm_q_vec_init; zeros(N_est,1)];
        else
            t_meas(k+1) = time_prop_utc_ekf(ii);
            ode_step = [t_meas(k),t_meas(k+1)];
            X_aug(:,1) = [X_ref(:,k);  stm_vec_init; stm_q_vec_init; zeros(N_est,1)];
            X1 = ode4(@rho_integ,ode_step,X_aug(:,1),parameters,Tsamp);        % run offline trajectory
            %                 [t,X1] = ode45(@(t,x)propagator_dmc(t, x, parameters), ode_step, X_aug(:,1), options);
            X_aug = X1(end,:)';
            
        end
        X_ref(:,k+1) = X_aug(1:N_st);
        stm_vec_full = X_aug(N_st+1:N_st+N_est^2);
        stm_rho = reshape(stm_vec_full,N_est,N_est);
        stm_f = eye(N_f);
        stm = [stm_rho, zeros(N_est, N_f);
            zeros(N_f, N_est), stm_f];
        
        stm_q_vec = X_aug(N_st+N_est^2+1:N_st+2*N_est^2);
        b_corr(:,k+1) = [X_aug(N_st+2*N_est^2+1:end);zeros(N_f,1)];
        if estimated_coeff.rho_DMC || parameters.empirical  && k > 0 && t_meas(k+1) - t_meas(k) < 100
            Q_rho = reshape(stm_q_vec,N_est,N_est);
            Q = [Q_rho, zeros(N_est, N_f);
                zeros(N_f, N_all)];
        else
            Q = zeros(N_est);
        end
        
        rho_est(ii) = rho_nom(ii)*(1+X_ref(1,k+1)+X_ref(3,k+1));
        y_pred(:,ii) =  -0.5*rho_est(ii)*Vi^2/(mass)*[Tx;Ty;Tz];
        
        Y_tor(:,ii) = y_meas(:,ii) - y_pred(:,ii);
        residual = Y_tor(:,ii);
        
        H_f = -0.5*rho_est(ii)*Vi^2/(mass)*[Xf_arg1(vec_est1)', zeros(1,N_est2), zeros(1,N_est3); zeros(1,N_est1), Xf_arg2(vec_est2)', zeros(1,N_est3);...
            zeros(1,N_est1),zeros(1,N_est2),Xf_arg3(vec_est3)'];
        
        H_rho = -0.5*Vi^2/(mass)*rho_nom(ii)*[Tx;Ty;Tz]*[1,0,1];
        H_tilde = [H_rho(:,vec_est), H_f];
        Omg = eye(N_all);
        [x_up, P_up, innov, P_pred, post_fit] = NLKF_cc(stm, H_tilde, Q, R_aug, Omg, residual, x_prior_fixed, P_up, vec_con);
        
        y_res_postfit(:,k+1) = post_fit;
        sigma_x(:,k+1) = sqrt(diag(P_up));
        
        x_up_rho = x_up(1:N_est);
        x_up_f = x_up(N_est+1:end);
        X_est_rho = X_ref(vec_est,k+1) + x_up_rho;
        X_ref(vec_est,k+1) = X_est_rho;
        X_est1(vec_est1) = X_est1(vec_est1) + x_up_f(1:N_est1);
        X_est2(vec_est2) = X_est2(vec_est2) + x_up_f(N_est1+1:N_est1+N_est2);
        X_est3(vec_est3) = X_est3(vec_est3) + x_up_f(N_est1+N_est2+1:N_est1+N_est2+N_est3);
        Xf_est = [X_est1;X_est2;X_est3];
        X_est(:,k+1) = [X_est_rho; X_est1(vec_est1); X_est2(vec_est2); X_est3(vec_est3)];
        
        
        
        
        Xest_mat(:,ii) = Xf_est;
        residual_mat(iter).residual(:,ii) = residual;
        
        %% Smoother storing
        smoother(k+1).stm = stm;
        smoother(k+1).Q = Q;
        smoother(k+1).P_pred = P_pred;
        smoother(k+1).P_up = P_up;
        smoother(k+1).Pxc_up = P_up;
        smoother(k+1).Pcc = P_pred;
        smoother(k+1).x_up = X_est(:,k+1);
        smoother(k+1).Gam = Omg;
        smoother(k+1).res = residual;
        smoother(k+1).H_tilde = H_tilde;
        smoother(k+1).bkk = b_corr(:,k+1);
        k = k+1;
    end
    %% Smoother
    
    kk = numel(t_meas);
    Xs_est(:,kk) = smoother(kk).x_up;
    Ps_up = smoother(kk).P_up;
    % error_xs(:,kk) = error_x(:,kk);
    ys_res_postfit(:,kk) = y_res_postfit(:,kk);
    meas_res_smooth(:,kk) = ys_res_postfit(:,kk);
    sigma_xs(:,kk) = sqrt(diag(Ps_up));
%     tr_Ps(kk) = Ps_up(1,1)+Ps_up(2,2)+Ps_up(3,3);
%     tr_Vs(kk) = Ps_up(4,4)+Ps_up(5,5)+Ps_up(6,6);
    
    Xs_est(:,kk) = X_est(:,kk);
    rhos_est(kk) = rho_est(kk);
    %     error_xs(1:6,kk) = error_x(1:6,kk);
    for kk = numel(t_meas)-1:-1:1
        
        phi = smoother(kk+1).stm;
        Pkk_up = smoother(kk).P_up;
        Pkk_pred = smoother(kk+1).P_pred;
        Gamkk = smoother(kk).Gam;
        Qk = smoother(kk).Q;
        xkk = smoother(kk).x_up;
        bkk = smoother(kk+1).bkk;
        P_low = chol(Pkk_pred, 'lower');
        Lambda = inv(P_low')*inv(P_low);
        Sk = Pkk_up*phi'*Lambda;
        
        
        Xs_est(:,kk) = xkk + Sk*(Xs_est(:,kk+1) -phi*xkk-bkk);
        Ps_up = Pkk_up + Sk*(Ps_up - Pkk_pred)*Sk';
        
        sigma_xs(:,kk) = sqrt(diag(Ps_up));
        
%         tr_Ps(kk) = Ps_up(1,1)+Ps_up(2,2)+Ps_up(3,3);
%         tr_Vs(kk) = Ps_up(4,4)+Ps_up(5,5)+Ps_up(6,6);
        
%         ys_res_postfit(:,kk) = y_meas(:,kk)- Xs_est(1:6,kk);
%         meas_res_smooth(:,kk) = ys_res_postfit(:,kk);
        X_ref(vec_est,kk) = Xs_est(1:N_est,kk);
        
        rhos_est(kk) = rho_nom(kk)*(1+X_ref(1,kk)+X_ref(3,kk));
    end

    x0_est = Xs_est(1:N_est,kk);
    Xs_f = Xs_est(N_est+1:end,kk);
    xcorr = [X_nom(vec_est) - Xs_est(1:N_est,1); X_nom_f - Xs_f]
    
    
    X_est1 = Xf_est(1:N_f1); X_est2 = Xf_est(N_f1+1:N_f1+N_f2); X_est3 = Xf_est(N_f1+N_f2+1:end);
    X_est1(vec_est1) = Xs_f(1:N_est1);
    X_est2(vec_est2) = Xs_f(N_est1+1:N_est1+N_est2);
    X_est3(vec_est3) = Xs_f(N_est1+N_est2+1:N_est1+N_est2+N_est3);
    Xf_est = [X_est1;X_est2;X_est3];
    x0_est_iter(:,iter) = xcorr;
    X_nom(vec_est) = x0_est;
    X_nom_f = Xs_f;
    rms_res_post(:,iter) = rms(y_res_postfit');
    rms_res_post_sm(:,iter) = rms(ys_res_postfit');
    % Cd calc
    if norm(xcorr) < 1e-3
        break
    end
    
end

err_final = [X_f1t(vec_est1); X_f2t(vec_est2); X_f3t(vec_est3)] - [X_est1(vec_est1); X_est2(vec_est2); X_est3(vec_est3)];

%% Plots
figure()
subplot(3,1,1)
plot(time_prop_utc/3600,residual_mat(1).residual(1,:), 'LineWidth',2)
hold on
plot(time_prop_utc/3600,residual_mat(iter).residual(1,:), 'LineWidth',2)
ylabel('X (Nm)')
title('Torque estimation errors')
legend('Initial','Estimation')
set(gca,'FontSize',16)
subplot(3,1,2)
plot(time_prop_utc/3600,residual_mat(1).residual(2,:), 'LineWidth',2)
hold on
plot(time_prop_utc/3600,residual_mat(iter).residual(2,:), 'LineWidth',2)
ylabel('Y (Nm)')
set(gca,'FontSize',16)
subplot(3,1,3)
plot(time_prop_utc/3600,residual_mat(1).residual(3,:), 'LineWidth',2)
hold on
plot(time_prop_utc/3600,residual_mat(iter).residual(3,:), 'LineWidth',2)
ylabel('Z (Nm)')
xlabel('Time (hours)')
set(gca,'FontSize',16)


figure()
plot(time_prop_utc/3600,rho, 'LineWidth',2)
hold on
plot(time_prop_utc/3600,rho_nom, 'LineWidth',2)
plot(time_prop_utc/3600,rhos_est, 'LineWidth',2)
ylabel('Density ($kg/m^3$)','Interpreter','latex')
xlabel('Time (hours)')
title('Density')
legend('True (HASDM)','Filter baseline (NRLMSISE-00)','Estimate')
grid on
set(gca,'FontSize',16)

%% Bias decorrelation
load 300km_torquesimulestJB08_iter1
figure()
plot(time_prop_utc/3600,rho,'k', 'LineWidth',2)
hold on
plot(time_prop_utc/3600,rho_nom, 'LineWidth',2)
plot(time_prop_utc/3600,rhos_est,'g', 'LineWidth',1)
hold on
load 300km_torquesimulestJB08_iter2
plot(time_prop_utc/3600,rhos_est, '--m', 'LineWidth',1)
ylabel('Density ($kg/m^3$)','Interpreter','latex')
xlabel('Time (hours)')
title('Density')
legend('True (HASDM)','Filter baseline (NRLMSISE-00)','Estimate (Iteration 1)', 'Estimate (Iteration 2)')
grid on
set(gca,'FontSize',16)