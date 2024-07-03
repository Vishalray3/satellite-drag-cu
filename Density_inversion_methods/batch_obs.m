%% batch estimator for observability study
tic
%% Batch filter algorithm
rms_residual(1) = 100;
rms_residual(2) = 1;
tol = 1e-2;
iter = 2;
x0_est = zeros(N_st,1);
P0_est = P_prior;
P_sqrt = sqrtm(P_prior);
N_iter = 3; %5;
step = time_prop_utc;
tym = step;
parameters.tym = tym;
for iter = 1:N_iter
    % while abs(rms_residual(iter) - rms_residual(iter-1))>tol
    iteration = iter-2
    rms_res = 0;
    P_up = P0_est;
    sigma_x(:,1) = sqrt(diag(P_up));
    X_nom = X_nom + x0_est;
    x_prior = x_prior - x0_est;
    P_low = chol(P_prior, 'lower');
    Lambda = inv(P_low')*inv(P_low); %diag([0.1,0.1,0.1,1e2,1e2,1e2,zeros(1,Order_o+1)]); %  % zeros(N_st,N_st); %                     % initial information matrix
    obs_mat = zeros(N_st,N_st);
    Hm = zeros(6,N_st);
    N =   Lambda*x_prior; % zeros(N_st,1); %
    k=0;
    X_aug(:,1) = [X_nom(:,1); stm_init(:)];
    X1 = ode4(@propagator_num,ode_step,X_aug(:,1),parameters,Tsamp);
    X_aug = X1';
    for i = 1:numel(time_prop)
        if isnan(time_prop(i))
            continue;
        else
            t_meas(k+1) = time_prop(i);
            X_ref(:,k+1) = X_aug(1:N_st,i);
            stm_vec = X_aug(N_st+1:N_st+N_st^2,i);
            stm = reshape(stm_vec,N_st,N_st);
%             obs_vec = X_aug(N_st+N_st^2+1:end,i);
%             obs_mat = reshape(obs_vec,N_st,N_st);
%             s_obs(:,k+1) = svd(obs_mat);
%             s_thresh(:,k+1) = N_st*max(s_obs(:,k+1))*2.220446049250313e-16;
%             STM(i).stm = stm;
            
            [y_pred, H_tilde] = meas_GPS( t_meas(k+1), X_ref(:,k+1), meas,N_st);
%             H_mat(i).H_tilde = H_tilde;
%             dely(:,k+1) = H_tilde(:,1:6)*stm(1:6, N_nom+1:N_st)*delp;
            H = H_tilde*stm;
            Hm = Hm+H;
            H_mat_cum(:,:,i) = Hm;
            H_mat(:,:,i) = H;
            residual = y_meas(:,i) - y_pred;
            y_res(:,k+1) = residual;
            Lambda = Lambda + H'/R_aug*H;
            obs_rank_crlb(k+1) = rank(Lambda);
            s_obs_crlb(:,k+1) = svd(Lambda);
            s_thresh_crlb(:,k+1) = N_st*max(s_obs_crlb(:,k+1))*2.220446049250313e-16;     
            cond_crlb(k+1) = max(s_obs_crlb(:,k+1))/min(s_obs_crlb(:,k+1));
            
            obs_mat = obs_mat+H'/R_aug*H;
            
            obs_rank(k+1) = rank(obs_mat);
            s_obs(:,k+1) = svd(obs_mat);
            s_thresh(:,k+1) = N_st*max(s_obs(:,k+1))*2.220446049250313e-16;
            cond(k+1) = max(s_obs(:,k+1))/min(s_obs(:,k+1));
%             obs_rank_con(k+1) = rank(obs_mat);
%             P_up = stm*P0_est*stm';
            Lambda_low = chol(Lambda,'lower');
            Lambda_inv = inv(Lambda_low')*inv(Lambda_low);
            P_up = Lambda_inv;
            sigma_x(:,k+1) = sqrt(diag(P_up));

%             P_up_nor = P_sqrt\P_up/P_sqrt;
%             P_up_nor = N_st/trace(P_up_nor)*P_up_nor;
%             [V,D] = eig(P_up_nor);
            err_x(:,k+1) = -X_true(:,k+1) + X_aug(1:6,k+1);
            err_stm(:,k+1) = stm(1:6,1:6)*del_X;                           % add errors to parameters and check how the errors change
            del_err(:,k+1) = abs(err_x(:,k+1) - err_stm(:,k+1));
            tr_P(k+1) = P_up(1,1)+P_up(2,2)+P_up(3,3);
            tr_V(k+1) = P_up(4,4)+P_up(5,5)+P_up(6,6);
            N = N + H'/R_aug*residual;
            rms_res = rms_res + residual'/R_aug*residual;
            k = k+1;
        end
    end
    iter = iter+1;
    rms_res_post(:,iter) = rms(y_res');
    Lambda_low = chol(Lambda,'lower');
    Lambda_inv = inv(Lambda_low')*inv(Lambda_low);
    X_est = X_ref;
    x0_est = Lambda_inv*N
    P0_est = Lambda_inv;
    xhat_iter(:,iter) = x0_est;
    rms_residual(iter) = rms_res;
    
end
%% Post-processing
sigma_mat = repmat(sigma_x(:,end),1,N_st);
Corr_mat = P0_est./(sigma_mat.*sigma_mat');

