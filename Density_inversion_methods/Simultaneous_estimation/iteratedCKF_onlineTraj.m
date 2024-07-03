%% Classical Kalman filter
fprintf("Running CKF")

%%
tr_P = NaN(1,numel(time_prop));
tr_V = NaN(1,numel(time_prop));
tr_Ps = NaN(1,numel(time_prop));
tr_Vs = NaN(1,numel(time_prop));
y_res_postfit = NaN(numel(meas),numel(time_prop));
y_res_prefit = NaN(numel(meas),numel(time_prop));
meas_res_smooth = NaN(numel(meas),numel(time_prop));

% X_est = X_nom;
N_iter = 6;

% str_est = [{'A20'},{'A40'},{'B10'}];
% ind_est = find(ismember(Xf_name_all,str_est));
% vec_est = [1:10,ind_est'+9];
vec_est = [1:8,11:numel(X_nom)];
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
P_prior = P_prior(vec_est,vec_est);
parameters.vec_est = vec_est;
X_est = zeros(N_est,numel(time_prop));
for iter=1:N_iter
    iteration = iter-1
    X_nom(vec_est) = X_nom(vec_est) + x0_est;
    x_prior = x_prior - x0_est;

    
    P_up = P_prior;
    sigma_x(:,1) = sqrt(diag(P_up));
    x_up = x_prior;

    for i = 1:numel(time_prop_utc)
        
        if isnan(time_prop_utc(i))
            continue;
        else
            if time_prop_utc(i) == 0
                k = 0;
                X_aug(:,1) = [X_nom(:,1);  stm_vec_init; stm_q_vec_init];
            else
                t_meas(k+1) = time_prop_utc(i);
                ode_step = [t_meas(k),t_meas(k+1)];
                X_aug(:,1) = [X_ref(:,k);  stm_vec_init; stm_q_vec_init];
                X1 = ode4(@propagator_dmc,ode_step,X_aug(:,1),parameters,Tsamp);        % run offline trajectory
                X_aug = X1(end,:)';
            end
            X_full(:,k+1) = X_aug;
            X_ref(:,k+1) = X_aug(1:N_st);
            stm_vec_full = X_aug(N_st+1:N_st+N_est^2);
            stm = reshape(stm_vec_full,N_est,N_est);
            stm_q_vec = X_aug(N_st+N_est^2+1:end);
            if ismember('rho_DMC', estimated_coeff) || ismember('cd_DMC', estimated_coeff)
                Q = reshape(stm_q_vec,N_est,N_est);
            else
                Q = zeros(N_est);
            end
            [y_pred(:,i), H_tilde] = meas_GPS( time_prop(i), X_ref(:,k+1), meas,N_est);

            residual = y_meas(:,i) - y_pred(:,i);
            %%%%%%%%%%%%%%%
            
            Omg = eye(N_est);
            [x_up,sig, P_up, innov, P_pred, post_fit] = NLKF_cc(stm, H_tilde, Q, R_aug, Omg, residual, x_up, P_up, vec_con);
            
            y_res_prefit(:,k+1) = innov;
            
            y_res_postfit(:,k+1) = post_fit;
            sigma_x(:,k+1) = sig;
%             sigma_c(:,k+1) = sig_c;
            X_est(:,k+1) = X_ref(vec_est,k+1) + x_up;
            
            error_x(1:6,k+1) = X_true(:,k+1) - X_est(1:6,k+1);
            tr_P(k+1) = P_up(1,1)+P_up(2,2)+P_up(3,3);
            tr_V(k+1) = P_up(4,4)+P_up(5,5)+P_up(6,6);
            trP(k+1) = tr_P(k+1);
            trV(k+1) = tr_V(k+1);
            %% Smoother storing
            smoother(k+1).stm = stm;
            smoother(k+1).Q = Q;
            smoother(k+1).P_pred = P_pred;
            smoother(k+1).P_up = P_up;
            smoother(k+1).Pxc_up = P_up;
            smoother(k+1).Pcc = P_pred;
            smoother(k+1).x_up = x_up;
            smoother(k+1).Gam = Omg;
            smoother(k+1).res = residual;
            smoother(k+1).H_tilde = H_tilde;
            k = k+1;
        end
    end
    rms_res_post(:,iter) = rms(y_res_postfit');
    %% Smoother
    
    kk = numel(t_meas);
    xs_up(:,kk) = smoother(kk).x_up;
    Ps_up = smoother(kk).P_up;
    % error_xs(:,kk) = error_x(:,kk);
    ys_res_postfit(:,kk) = y_res_postfit(:,kk);
    meas_res_smooth(:,kk) = ys_res_postfit(:,kk);
    sigma_xs(:,kk) = sqrt(diag(Ps_up));
    tr_Ps(kk) = Ps_up(1,1)+Ps_up(2,2)+Ps_up(3,3);
    tr_Vs(kk) = Ps_up(4,4)+Ps_up(5,5)+Ps_up(6,6);
    
    Xs_est(:,kk) = X_est(:,kk);
    error_xs(1:6,kk) = error_x(1:6,kk);
    for kk = numel(t_meas)-1:-1:1
        
        phi = smoother(kk+1).stm;
        Pkk_up = smoother(kk).P_up;
        Pkk_pred = smoother(kk+1).P_pred;
        Gamkk = smoother(kk).Gam;
        Qk = smoother(kk).Q;
        xkk = smoother(kk).x_up;
        P_low = chol(Pkk_pred, 'lower');
        Lambda = inv(P_low')*inv(P_low);
        Sk = Pkk_up*phi'*Lambda;
        

        xs_up(:,kk) = xkk + Sk*(xs_up(:,kk+1) -phi*xkk);
        Ps_up = Pkk_up + Sk*(Ps_up - Pkk_pred)*Sk';
        
        sigma_xs(:,kk) = sqrt(diag(Ps_up));
        Xs_est(:,kk) = X_ref(vec_est,kk) + xs_up(:,kk);
        error_xs(1:6,kk) = X_true(:,kk) - Xs_est(1:6,kk);
        
        tr_Ps(kk) = Ps_up(1,1)+Ps_up(2,2)+Ps_up(3,3);
        tr_Vs(kk) = Ps_up(4,4)+Ps_up(5,5)+Ps_up(6,6);
        
        ys_res_postfit(:,kk) = smoother(kk).res - smoother(kk).H_tilde*xs_up(:,kk);
        meas_res_smooth(:,kk) = ys_res_postfit(:,kk);
    end
    P_update(vec_est,vec_est) = Ps_up;
    x0_est = xs_up(:,1);
    %     x_corr(:,iter) = x0_est
    x0_est
    error_iter(iter).errorx = error_x;
    error_iter(iter).errorxs = error_xs;
%     if mod(iter,2) == 0 || iter == 1
%         X_iter(iter).x_est = X_est(7,:);
%         X_iter(iter).xs_est = Xs_est(7,:);
%         X_iter(iter).bias_s = Xs_est(9,:);
%         X_iter(iter).bias = X_est(9,:);
%     end
%     if iter < N_iter
%     clear xs_up Ps_up sigma_xs Xs_est Pxc_s 
%     end
    x0_est_iter(:,iter) = x0_est;
end
rms_res_post = rms(y_res_postfit');
rms_res_pre = rms(y_res_prefit');
rms_state = rms(error_x');

rms_pos = rms(rms_state(1:3));

rms_vel = rms(rms_state(4:6));

rms_res_sm = rms(ys_res_postfit');
rms_res_sm_all = rms(rms_res_sm/R_aug);
rms_res_all = rms(rms_res_post/R_aug);