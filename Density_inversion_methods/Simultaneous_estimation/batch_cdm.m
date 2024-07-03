%% batch estimator for observability study
clc
clearvars


y_meas(:,1) = [-172.653107;-880.207155;-6831.509509;-7.547016303;1.335101513;0.017315361];
% y_meas(:,2) = [-172.141938;-880.28864;-6831.516852;-7.547017665;1.335096758;0.01672375];
y_meas(:,2) = [-172.041938;-880.78864;-6831.516852;-7.547017665;1.335096758;0.01672375];

crr = 1360.761412841663;
ctr = -27980.65738302823;
ctt= 7035982.427392281;
cnr = 210.3805014183933;
cnt = 1826.946827301817; 
cnn = 331.0235183557189;
crdot_r = 31.24496325037386;
crdot_t = -7758.509120174653;
crdot_n = -2.021455208123943;
crdot_rdot = 8.556645538646126; 
ctdot_r = -1.403584331646452; 
ctdot_t = 10.13133280441347;
ctdot_n = -0.2357901627533037; 
ctdot_rdot = -0.01158842442023829;
ctdot_tdot = 0.001502579392117563;
cndot_r = 0.2335050550851493;
cndot_t = -2.265365070494777;
cndot_n = 0.1171002985748019;
cndot_rdot = 0.002942687318839407;
cndot_tdot = -0.0002453547281023429;
cndot_ndot = 0.0007498781740970351;
Rcov(:,:,1) = cov_mat_rtn(crr, ctr, ctt, cnr, cnt, cnn, crdot_r, crdot_t, crdot_n, crdot_rdot, ctdot_r, ctdot_t, ctdot_n, ctdot_rdot, ctdot_tdot, ...
    cndot_r, cndot_t, cndot_n, cndot_rdot, cndot_tdot, cndot_ndot)*1e-6;

R_rtn2eme = rtn2eme(y_meas(:,1));
y_meas(1:3,1) = R_rtn2eme'*y_meas(1:3,1); 
y_meas(4:6,1) = R_rtn2eme'*y_meas(4:6,1); 

crr = 1085.098354909838;
ctr = -10002.69759237779;
ctt= 1083219.802766544;
cnr = 147.5044262885784;
cnt = -37.44591344794553; 
cnn = 276.8748619932582;
crdot_r = 11.30786093655574;
crdot_t = -1192.172185295855;
crdot_n = 0.01092988358668493;
crdot_rdot = 1.313314348910292; 
ctdot_r = -1.154920605793044; 
ctdot_t = 6.258576173565409;
ctdot_n = -0.1635064470273625; 
ctdot_rdot = -0.00720096055584162;
ctdot_tdot = 0.001249400023241859;
cndot_r = 0.09747341476661531;
cndot_t = 0.5552305408743415;
cndot_n = 0.0729028348546945;
cndot_rdot = -0.0002478165135712488;
cndot_tdot = -0.0001048234101286761;
cndot_ndot = 0.0006455186541057117;
Rcov(:,:,2) = cov_mat_rtn(crr, ctr, ctt, cnr, cnt, cnn, crdot_r, crdot_t, crdot_n, crdot_rdot, ctdot_r, ctdot_t, ctdot_n, ctdot_rdot, ctdot_tdot, ...
    cndot_r, cndot_t, cndot_n, cndot_rdot, cndot_tdot, cndot_ndot)*1e-6;

R_rtn2eme = rtn2eme(y_meas(:,2));
y_meas(1:3,2) = R_rtn2eme'*y_meas(1:3,2); 
y_meas(4:6,2) = R_rtn2eme'*y_meas(4:6,2); 

time_prop = [0,0];


P_prior = diag(diag(Rcov(:,:,1)))*1e10;
x_prior = zeros(6,1);
%% Batch filter algorithm
rms_residual(1) = 100;
rms_residual(2) = 1;
tol = 1e-2;
iter = 2;


N_st = 3;
vec = 1:N_st;
P_prior = P_prior(vec,vec);
x_prior = x_prior(vec);
Rcov = Rcov(vec,vec,:);
y_meas = y_meas(vec,:);

x0_est = zeros(N_st,1);
P0_est = P_prior;
P_sqrt = sqrtm(P_prior);
N_iter = 3; %5;
stm_init = eye(N_st);
X_nom = y_meas(:,1);

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
    X_aug(:,2) = X_aug(:,1);
    for i = 1:numel(time_prop)
        if isnan(time_prop(i))
            continue;
        else
            t_meas(k+1) = time_prop(i);
            X_ref(:,k+1) = X_aug(1:N_st,i);
            stm_vec = X_aug(N_st+1:N_st+N_st^2,i);
            stm = reshape(stm_vec,N_st,N_st);

            y_pred = X_ref(vec,k+1);
            H_tilde = eye(N_st);
            
            H = H_tilde*stm;
            H_mat(:,:,i) = H;
            residual = y_meas(:,i) - y_pred;
            y_res(:,k+1) = residual;
            
            R_aug = Rcov(:,:,i);
            
%             R_aug = [Rcov(1:3,1:3,i), zeros(3,3);
%                     zeros(3,3), Rcov(4:6,4:6,i)];

%             V = chol(R_aug, 'lower');
%             H = V\H;
%             residual = V\residual;
%             R_aug = eye(6);

%             R_aug = diag(diag(Rcov(:,:,i)));
            
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

p = 1-10^(-1/2);

mu_est_ri = [0,0];
Sigma_est_ri = P_up(1:2,1:2);
mu_est_ci = [0,0];
Sigma_est_ci = P_up(2:3,2:3);

mu1_ri = y_meas(1:2,1) - X_nom(1:2);
Sigma1_ri = Rcov(1:2,1:2,1);
mu1_ci = y_meas(2:3,1) - X_nom(2:3);
Sigma1_ci = Rcov(2:3,2:3,1);

mu2_ri = y_meas(1:2,2) - X_nom(1:2);
Sigma2_ri = Rcov(1:2,1:2,2);
mu2_ci = y_meas(2:3,2) - X_nom(2:3);
Sigma2_ci = Rcov(2:3,2:3,2);

% mu_est_ri = [0,0];
% Sigma_est_ri = P_up(4:5,4:5);
% mu_est_ci = [0,0];
% Sigma_est_ci = P_up(5:6,5:6);
% 
% mu1_ri = y_meas(4:5,1) - X_nom(4:5);
% Sigma1_ri = Rcov(4:5,4:5,1);
% mu1_ci = y_meas(5:6,1) - X_nom(5:6);
% Sigma1_ci = Rcov(5:6,5:6,1);
% 
% mu2_ri = y_meas(4:5,2) - X_nom(4:5);
% Sigma2_ri = Rcov(4:5,4:5,2);
% mu2_ci = y_meas(5:6,2) - X_nom(5:6);
% Sigma2_ci = Rcov(5:6,5:6,2);
%%
figure()
subplot(2,1,1)
a = plotErrorEllipse(mu_est_ri, Sigma_est_ri, p);
plot(a(2, :) + mu_est_ri(2), a(1, :) + mu_est_ri(1),'k', 'LineWidth',1);
hold on
a = plotErrorEllipse(mu1_ri, Sigma1_ri, p);
plot(a(2, :) + mu1_ri(2), a(1, :) + mu1_ri(1),'b', 'LineWidth',1);
a = plotErrorEllipse(mu2_ri, Sigma2_ri, p);
plot(a(2, :) + mu2_ri(2), a(1, :) + mu2_ri(1),'r', 'LineWidth',1);
plot(mu_est_ri(2), mu_est_ri(1),'ko', 'MarkerFaceColor', 'k')
plot(mu1_ri(2), mu1_ri(1),'bo', 'MarkerFaceColor', 'b')
plot(mu2_ri(2), mu2_ri(1),'ro', 'MarkerFaceColor', 'r')


grid on
set(gca, 'FontSize', 14)
xlabel('In-track (km)')
ylabel('Radial (km)')
legend('Estimate', 'Observed','Observed')
title('RIC position w.r.t final estimate')

subplot(2,1,2)
a = plotErrorEllipse(mu_est_ci, Sigma_est_ci, p);
plot(a(1, :) + mu_est_ci(1), a(2, :) + mu_est_ci(2),'k', 'LineWidth',1);
hold on
plot(mu_est_ci(1), mu_est_ci(2),'ko', 'MarkerFaceColor', 'k')
a = plotErrorEllipse(mu1_ci, Sigma1_ci, p);
plot(a(1, :) + mu1_ci(1), a(2, :) + mu1_ci(2),'b', 'LineWidth',1);
plot(mu1_ci(1), mu1_ci(2),'bo', 'MarkerFaceColor', 'b')
a = plotErrorEllipse(mu2_ci, Sigma2_ci, p);
plot(a(1, :) + mu2_ci(1), a(2, :) + mu2_ci(2),'r', 'LineWidth',1);
plot(mu2_ci(1), mu2_ci(2),'ro', 'MarkerFaceColor', 'r')
grid on
set(gca, 'FontSize', 14)
xlabel('In-track (km)')
ylabel('Cross-track (km)')
%% plots
function a = plotErrorEllipse(mu, Sigma, p)

    s = -2 * log(1 - p);

    [V, D] = eig(Sigma * s);

    t = linspace(0, 2 * pi);
    a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
end

function Rcov = cov_mat_rtn(crr, ctr, ctt, cnr, cnt, cnn, crdot_r, crdot_t, crdot_n, crdot_rdot, ctdot_r, ctdot_t, ctdot_n, ctdot_rdot, ctdot_tdot, ...
    cndot_r, cndot_t, cndot_n, cndot_rdot, cndot_tdot, cndot_ndot)

ctr = ctr;
cnr = 0.0005*cnr;
cnt = 0.005*cnt;
% cnn = 1e10*cnn;
% crr = 1e2*crr;
% ctt = 1e10*ctt;
Rcov = 0.0005*[crr,     ctr,     cnr,     crdot_r,    ctdot_r,    cndot_r; 
        ctr,     ctt,     cnt,     crdot_t,    ctdot_t,    cndot_t;
        cnr,     cnt,     cnn,     crdot_n,    ctdot_n,    cndot_n;
        crdot_r, crdot_t, crdot_n, crdot_rdot, ctdot_rdot, cndot_rdot;
        ctdot_r, ctdot_t, ctdot_n, ctdot_rdot, ctdot_tdot, cndot_tdot;
        cndot_r, cndot_t, cndot_n, cndot_rdot, cndot_tdot, cndot_ndot];
end

function R_rtn2eme = rtn2eme(X)
R = X(1:3)/norm(X(1:3));
N = cross(R, X(4:6)/norm(X(4:6)));
N = N/norm(N);
T = cross(N,R);
T = T/norm(T);

R_rtn2eme = [R,T,N];
end