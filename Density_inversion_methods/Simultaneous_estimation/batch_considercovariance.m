%% batch estimator for consider covariance analysis
tic
%% Define parameters specific to consider covariance
  %cube bodf [7,15,8, 23, 31, 16]; % gpm bodf[1:7,9,11,13,28,30,68,92];    % states to be estimated
vec_con = [];%cube bodf []; % gpm bodf[8,10,12,14:27,29,31:67,69:91,93:N_st]; % states to be considered
N_est = numel(vec_est);
N_con = numel(vec_con);
parameters.vec_est = vec_est;
%% Batch filter algorithm
rms_residual(1) = 100;
rms_residual(2) = 1;
tol = 1e-2;
x0_est = zeros(N_est,1);
x_prior = x_prior(vec_est);
Pxx = P_prior(vec_est,vec_est);
N_iter = 1; %5;
step = time_prop_utc;
tym = step;
parameters.tym = tym;
% R_aug = R_aug(vec_meas,vec_meas);
% Initial covariances
Px_bar = Pxx;
Pcc_bar = P_prior(vec_con,vec_con);
Pxc_bar = zeros(N_est,N_con);
% stm_init = eye(N_est);
c_bar = X_nom(vec_con); %Af_err_init(vec_con-7)';%Xf_true(vec_con-7)';%zeros(numel(vec_con),1); %3*Xf_std(vec_con-7)';

% X_nom(8:end) = 0;
for iter = 1:N_iter
    % while abs(rms_residual(iter) - rms_residual(iter-1))>tol
    iteration = iter-2
    rms_res = 0;
    P_up = Pxx;
    sigma_x(:,1) = sqrt(diag(P_up));
    X_nom(vec_est) = X_nom(vec_est) + x0_est;
    x_prior = x_prior - x0_est;
    Px_low = chol(Px_bar, 'lower');             % covariance of the state
    Px_inv = inv(Px_low')*inv(Px_low);
    Pcc_low = chol(Pcc_bar, 'lower');           % covariance of consider parameters
    Pcc_inv = inv(Pcc_low')*inv(Pcc_low);
    Mxx_low = chol(Px_bar - Pxc_bar*Pcc_inv*Pxc_bar', 'lower');
    Mxx = inv(Mxx_low')*inv(Mxx_low);
    Mxc = -Mxx*Pxc_bar*Pcc_inv;
    Mcc_low = chol(Pcc_bar - Pxc_bar'*Px_inv*Pxc_bar,'lower');
    Mcc = inv(Mcc_low')*inv(Mcc_low);
    
    Nx =  Px_inv*x_prior; % zeros(N_st,1); %
    k=0;
    X_aug(:,1) = [X_nom(:,1); stm_init(:)];
    
%     [t,X1] = ode45(@(t,x)propagator_dmc(t, x, parameters), step, X_aug(:,1), options);
%     X_aug = X1';
    X1 = ode4(@propagator_dmc,ode_step,X_aug(:,1),parameters,Tsamp);
%     X_aug = interp1(ode_step, X1,sod_sec_sp3 - sod_sec_sp3(1),'spline');
%     X_aug = X_aug';
    X_aug = X1';
     
    for i = 1:numel(ode_step) %numel(sod_sec_sp3)
        if isnan(time_prop(i))
            continue;
        else
            t_meas(k+1) = time_prop(i);
            X_ref(:,k+1) = X_aug(1:N_st,i);
            stm_vec = X_aug(N_st+1:end,i);
            stm = reshape(stm_vec,N_st,N_st);
            %             STM(i).stm = stm;
            
            y_pred = X_ref(vec_meas,k+1);        % measurements are eci states
            H_tilde = [eye(6), zeros(6,N_st-6)];
            H_tilde = H_tilde(vec_meas,:);            
%             [y_pred, H_tilde] = meas_GPS( t_meas(k+1), X_ref(:,k+1), vec_meas,N_st);

            
            H = H_tilde*stm;
            Hx = H(:,vec_est);
            Hc = H(:,vec_con);
            Mxx = Mxx + Hx'/R_aug*Hx;
            Mxc = Mxc + Hx'/R_aug*Hc;
            Mcc = Mcc + Hc'/R_aug*Hc;
            
            
            residual = y_meas(:,i) - y_pred;
            y_res(:,k+1) = residual;
            Nx = Nx + Hx'/R_aug*residual;
            
            
            Mxx_low = chol(Mxx,'lower');
            Mxx_inv = inv(Mxx_low')*inv(Mxx_low);
            Px = Mxx_inv;                         % covariance from batch estimator that assumes no error, covariances in considered parameters
            Sxc = -Px*Mxc;
            P_up = Px + Sxc*Pcc_bar*Sxc';         % covariance that considers errors in considered parameters
            
            sigma_x(:,k+1) = sqrt(diag(Px));
            sigma_xx(:,k+1) = sqrt(diag(P_up));
            tr_P(k+1) = P_up(1,1)+P_up(2,2)+P_up(3,3);
            tr_V(k+1) = P_up(4,4)+P_up(5,5)+P_up(6,6);
            rms_res = rms_res + residual'/R_aug*residual;
            k = k+1;
        end
    end
    iter = iter+1;
    rms_res_post(:,iter) = rms(y_res');
    Mxx_low = chol(Mxx,'lower');
    Mxx_inv = inv(Mxx_low')*inv(Mxx_low);
    Px = Mxx_inv;
    Sxc = -Px*Mxc;
    X_est = X_ref;
    if ~isempty(c_bar)
        x0_est = Mxx_inv*Nx - Mxx_inv*Mxc*c_bar;
        Pxx = Px + Sxc*Pcc_bar*Sxc';
    else
        x0_est = Mxx_inv*Nx
        Pxx = Px;
    end
    
    Pxc = Sxc*Pcc_bar;
    xhat_iter(:,iter) = x0_est;
    rms_residual(iter) = rms_res;
    
end