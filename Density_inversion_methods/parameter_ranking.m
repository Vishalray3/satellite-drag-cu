%% Parameter ranking by successive orthogonalization
% Ref: Methods for parameter ranking in nonlinear mechanistic models, Lund
% et al., 2005
clc
clearvars

addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Observability study/Results/BFF')
load('bodf_cube_all')
par_ind = 0:4:30; %81 1:102;

W1 = sqrt(R_aug); % diagonal matrix normalizing states
% W2 = diag(abs([Cd,X_f])); %diag(abs([del_X',Af0_true,X_f])); % diagonal matrix normalizing parameters

% stm_mat = [];
% for jj = 1:numel(time_pred_utc)
%     stm_vec = X_aug(N_st+1:N_st+N_st^2,jj);
%     stm = reshape(stm_vec,N_st,N_st);
%     stm_mat = [stm_mat;W1*stm(1:3,7:end)*W2];
% %     stm_sen = W1*stm(1:6,7:end)*W2;
% if jj > 10
%     [Q_sen,R_sen,P_sen] = qr(stm_mat,0);
%     rank_par(:,jj-10) = par_ind(P_sen);
%     norm_par(:,jj-10) = diag(R_sen).^2;
% end
% end

% W2 = eye(N_st);
W2 = diag(abs([Cd_std,Xf_std]));
H_sen_mat = [];
N_par = 6;

for jj = 1:numel(time_pred_utc)
    H_sen = H_mat(:,:,jj);
    H_sen_mat = [H_sen_mat;W1\H_sen(:,1+N_par:end)*W2];      % W2(N_par+1:end,N_par+1:end)
    
    %     if jj > 10
    %         [Q_sen,R_sen,P_sen] = qr(H_sen_mat,0);
    %         rank_par(:,jj-10) = par_ind(P_sen);
    %         norm_par(:,jj-10) = abs(diag(R_sen));
    %     end
end

[Q_sen,R_sen,P_sen] = qr(H_sen_mat,0);
rank_par_full = par_ind(P_sen);
norm_par_full = abs(diag(R_sen));
% Hsen_low = chol(H_sen_mat'*H_sen_mat,'lower');
% Hsen_inv = inv(Hsen_low')*inv(Hsen_low);
Hsen_inv = inv(R_sen'*R_sen);
sig_Hsen1 = diag(Hsen_inv);
Hsen_inv = inv(H_sen_mat'*H_sen_mat);
sig_Hsen2 = diag(Hsen_inv);
%% Individual variance contribution
% % diagonalize R_sen
% R_D = diag(diag(R_sen));
% R_bar = R_D\R_sen;
% U_bar = inv(R_bar);
% Hsen_inv = U_bar/(R_D^2)*U_bar';
% sig_Hsen3 = diag(Hsen_inv);
% 
% Unorm = (U_bar'*U_bar)/R_D^2;
% indiv_var = diag(Unorm);
% 
% Ubar2 = inv(chol(H_sen_mat'*H_sen_mat,'upper'));
% indiv_var2 = diag(Ubar2'*Ubar2);
% 
% %%
% vec = [1:10];
% H_low = chol(H_sen_mat(:,vec)'*H_sen_mat(:,vec)  ,'lower'); % +inv(P_prior(vec,vec))
% H_inv = inv(H_low')*inv(H_low);
% Cov_mat = H_inv;
% var_H = diag(Cov_mat);
% var_mat = repmat(sqrt(var_H),1,numel(vec));
% Corr_rank = Cov_mat./(var_mat.*var_mat');
%% Figure
% for ii = 1:numel(par_ind)
%     str_x{ii} = strcat('$\overline{A}_{',num2str(par_ind(ii)),'}$');
% end
% gpm
%%%%%%
% Xf_all = [{'A00'},Xf_name];
Xf_rank = Xf_name_all(rank_par_full);
% Xf_rank = regexprep(Xf_rank, '.$', '', 'lineanchors');
% str_x = [{'$\overline{A}_{0}$'},{'$\overline{A}_{2}$'},{'$\overline{B}_{1}$'},{'$\overline{A}_{4}$'},{'$\overline{B}_{3}$'},...
%     {'$\overline{A}_{6}$'},{'$\overline{A}_{8}$'},{'$\overline{A}_{1}$'},{'$\overline{A}_{10}$'},{'$\overline{A}_{12}$'},{'$\overline{A}_{14}$'},...
%     {'$\overline{A}_{16}$'},{'$\overline{A}_{18}$'},{'$\overline{B}_{5}$'},{'$\overline{A}_{3}$'},{'$\overline{A}_{20}$'},{'$\overline{A}_{22}$'},...
%     {'$\overline{A}_{24}$'},{'$\overline{B}_{7}$'},{'$\overline{A}_{26}$'}]
%gpm bodf
% str_x = [{'\boldmath$A_{0,0}$'},{'\boldmath$A_{2,0}$'},{'\boldmath$B_{1,0}$'},{'\boldmath$A_{4,0}$'},{'\boldmath$A_{0,1}$'},...
%     {'\boldmath$B_{3,0}$'},{'\boldmath$A_{6,0}$'},{'\boldmath$A_{2,1}$'},{'\boldmath$A_{8,0}$'},{'\boldmath$A_{10,0}$'},{'\boldmath$A_{12,0}$'},...
%     {'\boldmath$B_{1,1}$'},{'\boldmath$A_{14,0}$'},{'\boldmath$A_{16,0}$'},{'\boldmath$A_{2,3}$'},{'\boldmath$A_{18,0}$'},{'\boldmath$A_{20,0}$'},...
%     {'\boldmath$B_{3,1}$'},{'\boldmath$A_{22,0}$'},{'\boldmath$A_{24,0}$'},{'\boldmath$B_{5,0}$'},{'\boldmath$A_{6,1}$'},{'\boldmath$A_{8,1}$'}];

% cube bodf
% str_x = [{'$\overline{A}_{0,0}$'},{'$\overline{A}_{0,1}$'},{'$\overline{A}_{4,0}$'},{'$\overline{A}_{8,0}$'},{'$\overline{A}_{0,2}$'},...
%              {'$\overline{A}_{12,0}$'}, {'$\overline{A}_{0,3}$'}, {'$\overline{A}_{16,0}$'}, {'$\overline{A}_{20,0}$'}, {'$\overline{A}_{4,1}$'}];

% str_x = [{'$\overline{A}_{0}$'},{'$\overline{A}_{1}$'},{'$\overline{A}_{2}$'},{'$\overline{A}_{3}$'},{'$\overline{A}_{4}$'},...
%              {'$\overline{A}_{5}$'}, {'$\overline{A}_{6}$'}, {'$\overline{A}_{7}$'}, {'$\overline{A}_{8}$'}, {'$\overline{A}_{9}$'}, {'$\overline{A}_{10}$'}];
str_x = [{'$\overline{A}_{0}$'},{'$\overline{A}_{4}$'},{'$\overline{A}_{8}$'},{'$\overline{A}_{12}$'},{'$\overline{A}_{16}$'},...
             {'$\overline{A}_{20}$'}, {'$\overline{A}_{24}$'}, {'$\overline{A}_{28}$'}];

%%%%%%%%
figure()
h = bar(norm_par_full(1:end))
set(gca,'yscale','log')
set(gca,'FontSize',18)
ylabel('Sensitivity norm')
title('Orthogonal sensitivities')
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTickLabel', str_x);
yticks([1e-4, 1e-2, 1e0, 1e2, 1e4,1e6])
% yticks([1e-10, 1e-5, 1e0, 1e5])
ylim([1-10 1e7])
% xticks([1:23])
% xtickangle(45)
% ylim([1e-10 1e7])
grid on
ax = ancestor(h, 'axes');
Xrule = ax.XAxis;
% Xrule.FontSize = 12;
