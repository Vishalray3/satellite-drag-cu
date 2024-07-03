function [y_res_postfit, ys_respostfit, X_est, Xs_est, sigma_x, sigma_xs] = ekf_dmc_func(time_prop_utc_ekf,vec_est,P_prior,R_aug,vec_meas,X_nom,y_meas,...
                                                                            parameters,N_st,estimated_coeff,Tsamp)
%#codegen
tr_P = NaN(1,numel(time_prop_utc_ekf));
tr_V = NaN(1,numel(time_prop_utc_ekf));
tr_Ps = NaN(1,numel(time_prop_utc_ekf));
tr_Vs = NaN(1,numel(time_prop_utc_ekf));
y_res_postfit = NaN(numel(vec_meas),numel(time_prop_utc_ekf));
y_res_prefit = NaN(numel(vec_meas),numel(time_prop_utc_ekf));
meas_res_smooth = NaN(numel(vec_meas),numel(time_prop_utc_ekf));

% X_est = X_nom;
N_iter = 5;

% str_est = [{'A10'},{'A20'},{'A30'},{'A50'},{'A70'},{'B10'}];
% str_est = [{'A01'},{'A02'},{'C01'}]; %,{'A02'},{'C01'}];
% ind_est = find(ismember(Xf_name_all,str_est));
% vec_est = [1:9,ind_est'+9];
% vec_est = [1:6];%[1:9,11:numel(X_nom)];
vec_con = [];
N_est = numel(vec_est);
N_con = numel(vec_con);
x0_est = zeros(N_est,1);
P_update = P_prior;

stm_init = eye(N_est);
stm_q_init = zeros(N_est);
stm_vec_init = stm_init(:);
stm_q_vec_init = stm_q_init(:);
x_prior = zeros(N_est,1);
x_prior_fixed = x_prior;
P_prior = P_prior(vec_est,vec_est);
parameters.vec_est = vec_est;

X_est = zeros(N_est,numel(time_prop_utc_ekf));
R_aug = R_aug(vec_meas,vec_meas);
X_aug = zeros(numel(X_nom)+2*N_est^2+N_est,1);
t_meas = zeros(1,numel(time_prop_utc_ekf));
X_ref = zeros(N_st,numel(t_meas));
X_full = zeros(numel(X_aug),numel(t_meas));
sigma_x = zeros(N_est,numel(t_meas));
b_corr = zeros(N_est,numel(t_meas));
y_pred = zeros(numel(vec_meas),numel(t_meas));
smoother = struct();
ys_res_postfit = NaN(numel(vec_meas),numel(time_prop_utc_ekf));
Xs_est = NaN(N_est,numel(time_prop_utc_ekf));
sigma_xs = zeros(N_st,numel(t_meas));
x0_est_iter = zeros(N_est,N_iter);

for iter=1:N_iter
    iteration = iter-1
    
    P_up = P_prior;
    sigma_x(:,1) = sqrt(diag(P_up));
    
    for i = 1:numel(time_prop_utc_ekf)
        
        if isnan(time_prop_utc_ekf(i))
            continue;
        else
            if time_prop_utc_ekf(i) == 0
                k = 0;
                X_aug(:,1) = [X_nom(:,1);  stm_vec_init; stm_q_vec_init; zeros(N_est,1)];
            else
                t_meas(k+1) = time_prop_utc_ekf(i);
                if t_meas(k+1) - t_meas(k) <= 10
                    ode_step = [t_meas(k),t_meas(k+1)];
                    X_aug(:,1) = [X_ref(:,k);  stm_vec_init; stm_q_vec_init; zeros(N_est,1)];
                    X1 = ode4(@propagator_dmc,ode_step,X_aug(:,1),parameters,Tsamp);        % run offline trajectory
                    %                 [t,X1] = ode45(@(t,x)propagator_dmc(t, x, parameters), ode_step, X_aug(:,1), options);
                    X_aug = X1(end,:)';
                else
                    ode_step = [t_meas(k):10:t_meas(k+1)+10];
%                     X_ref(8:10,k) = 0;
                    X_aug(:,1) = [X_ref(:,k);  stm_vec_init; stm_q_vec_init; zeros(N_est,1)];
                    X1 = ode4(@propagator_dmc,ode_step,X_aug(:,1),parameters,Tsamp);        % run offline trajectory
                    %                 [t,X1] = ode45(@(t,x)propagator_dmc(t, x, parameters), ode_step, X_aug(:,1), options);
                    X_aug = interp1(ode_step, X1,t_meas(k+1),'spline');
                    X_aug = X_aug';
                end
                
            end
            X_full(:,k+1) = X_aug;
            X_ref(:,k+1) = X_aug(1:N_st);
            stm_vec_full = X_aug(N_st+1:N_st+N_est^2);
            stm = reshape(stm_vec_full,N_est,N_est);
            stm_q_vec = X_aug(N_st+N_est^2+1:N_st+2*N_est^2);
            b_corr(:,k+1) = X_aug(N_st+2*N_est^2+1:end);
            if estimated_coeff.rho_DMC  && k > 0 && t_meas(k+1) - t_meas(k) < 100
                Q = reshape(stm_q_vec,N_est,N_est);
            else
                Q = zeros(N_est);
            end
            y_pred(:,i) = X_ref(vec_meas,k+1);        % measurements are eci states
            H_tilde = [eye(6), zeros(6,N_est-6)];
            H_tilde = H_tilde(vec_meas,:);
            %             [y_pred(:,i), H_tilde] = meas_GPS( time_prop(i), X_ref(:,k+1), vec_meas,N_est);
            residual = y_meas(:,i) - y_pred(:,i);
            %%%%%%%%%%%%%%%
            
            Omg = eye(N_est);
            [x_up, P_up, innov, P_pred, post_fit] = NLKF_cc(stm, H_tilde, Q, R_aug, Omg, residual, x_prior_fixed, P_up, vec_con);
            
            y_res_prefit(:,k+1) = innov;
            
            y_res_postfit(:,k+1) = post_fit;
            sigma_x(:,k+1) = sqrt(diag(P_up));
%             crosscov_rho(k+1) = P_up(Nrho-2,Nrho);
            X_est(:,k+1) = X_ref(vec_est,k+1) + x_up;
            X_ref(vec_est,k+1) = X_est(:,k+1);
            %             error_x(1:6,k+1) = X_true(:,k+1) - X_est(1:6,k+1);
            tr_P(k+1) = P_up(1,1)+P_up(2,2)+P_up(3,3);
            tr_V(k+1) = P_up(4,4)+P_up(5,5)+P_up(6,6);

            
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
    end
    %% Smoother
    
    kk = numel(t_meas);
    Xs_est(:,kk) = smoother(kk).x_up;
    Ps_up = smoother(kk).P_up;
    % error_xs(:,kk) = error_x(:,kk);
    ys_res_postfit(:,kk) = y_res_postfit(:,kk);
    meas_res_smooth(:,kk) = ys_res_postfit(:,kk);
    sigma_xs(:,kk) = sqrt(diag(Ps_up));
    tr_Ps(kk) = Ps_up(1,1)+Ps_up(2,2)+Ps_up(3,3);
    tr_Vs(kk) = Ps_up(4,4)+Ps_up(5,5)+Ps_up(6,6);
    
    Xs_est(:,kk) = X_est(:,kk);
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
%         crosscov_rhos(k+1) = Ps_up(Nrho-2,Nrho);
        %         error_xs(1:6,kk) = X_true(:,kk) - Xs_est(1:6,kk);
        
        tr_Ps(kk) = Ps_up(1,1)+Ps_up(2,2)+Ps_up(3,3);
        tr_Vs(kk) = Ps_up(4,4)+Ps_up(5,5)+Ps_up(6,6);
        
        ys_res_postfit(:,kk) = y_meas(:,kk)- Xs_est(1:6,kk);
        meas_res_smooth(:,kk) = ys_res_postfit(:,kk);
        X_ref(vec_est,kk) = Xs_est(:,kk);
        
        
    end
    %     [Af_calc,Bf_calc] = gsim_estimation_func(Order_b,Order_o,Xf_name_all,X_ref(:,kk),Ps_up,parameters,time_prop(kk),oxa,Af_total,Bf_total);
    %     Cds_bias(kk) = Af_calc(1);
    %     parameters.Cd_est = Cds_bias(1);
    P_update(vec_est,vec_est) = Ps_up;
    x0_est = Xs_est(:,1);
    %     x_corr(:,iter) = x0_est
    xcorr = X_nom(vec_est) - Xs_est(:,1)
    %     error_iter(iter).errorx = error_x;
    %     error_iter(iter).errorxs = error_xs;
    %     if mod(iter,2) == 0 || iter == 1
    %         X_iter(iter).x_est = X_est(7,:);
    %         X_iter(iter).xs_est = Xs_est(7,:);
    %         X_iter(iter).bias_s = Xs_est(9,:);
    %         X_iter(iter).bias = X_est(9,:);
    %     end
    %     if iter < N_iter
    %     clear xs_up Ps_up sigma_xs Xs_est Pxc_s
    %     end
    x0_est_iter(:,iter) = xcorr;
    X_nom(vec_est) = x0_est;
%     rms_res_post(:,iter) = rms(y_res_postfit');
%     rms_res_post_sm(:,iter) = rms(ys_res_postfit');
    % Cd calc
    if norm(xcorr) < 1e-3
        break 
    end
end
clear smoother
% rms_res_post = rms(y_res_postfit');
% rms_res_pre = rms(y_res_prefit');
% % rms_state = rms(error_x');
% %
% % rms_pos = rms(rms_state(1:3));
% %
% % rms_vel = rms(rms_state(4:6));
%
% rms_res_sm = rms(ys_res_postfit');
% rms_res_sm_all = rms(rms_res_sm/R_aug);
% rms_res_all = rms(rms_res_post/R_aug);