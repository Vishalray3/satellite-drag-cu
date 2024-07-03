load('case10_gsimEst','Af_pred','Bf_pred')
load('case10_iter1','Xs_est','sigma_xs')
% load('case8_truth','X_true_aug')
Af_total_mat(:,:,K_ind) = Af_pred;
Bf_total_mat(:,:,K_ind) = Bf_pred;
del_r = sigma_xs(1:3,1)*10; %(Xs_est(1:3,1) - X_true_aug(1:3,1))*5; %[10;10;10];     % initial standard deviations
del_v = sigma_xs(4:6,1)*10; %(Xs_est(4:6,1) - X_true_aug(4:6,1))*5; %[0.1;0.1;0.1];

%% Initialize state

% X_init = [r_eci_true(:,1);v_eci_true(:,1)] + [del_r;del_v];
% 
% X_ref_st(:,1) = X_init;
X_init = [Xs_est(1:3,1);Xs_est(4:6,1)];%[reci(:,1);veci(:,1)];

X_ref_st(:,1) = X_init;
%% DMC model parameters
% X_s = Xs_est(7:9,1); 
% Xs_std = sigma_xs(7:9,1)'*10;
X_s = [0;0;Xs_est(8,1)]; 
Xs_std = [1;1;sigma_xs(8,1)']*10;
Cd_nom = Af_total_mat(1,1,K_ind);
clear X_true_aug Xs_est sigma_xs