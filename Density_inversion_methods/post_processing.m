%% post-processing
% addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Observability study/Results/OFF/Exponential')
% load('Case1_sphere')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Observability study/Results/BODF')
load('bodf_cube_all')
del_dash = area/mass;
A_effect(:,1) = zeros(11,1);
betaH = 1/H_scale;
N_orb = floor(time_pred_utc(end)/T_orb);
all_coeff = [Af0_true,Xf_true];
% for kk = 1:numel(all_coeff)     % Af_true for OFF
%     Af_coeff = [zeros(1,kk-1),all_coeff(kk)];     % Af_true for OFF
%     for jj = 1:numel(time_pred_utc) %N_orb
%         ind = jj; %ceil(jj*T_orb/del_T);
%         coe = rv2coe_E(X_est(1:3,ind),X_est(4:6,ind),mu_e);
%         a_est(jj) = coe(1);
%         e_est(jj) = coe(2);
%         ae_est(jj) = a_est(jj)*e_est(jj);
% %         [del_a,~] = analytical_change_off(del_dash,a_est(jj),e_est(jj),rho0,betaH,ae_est(jj),Af_coeff,kk-1,'low');
% %         [del_a,~] = analytical_change_bff_iner(del_dash,a_est(jj),e_est(jj),rho0,betaH,ae_est(jj),Af_coeff,kk-1,'low');
% %         A_effect(kk,jj+1) = A_effect(kk,jj) + del_a;
%
%         x_sen(kk,jj) = X_aug(N_st+(5+kk)*N_st+1,ind);
%         y_sen(kk,jj) = X_aug(N_st+(5+kk)*N_st+2,ind);
%         z_sen(kk,jj) = X_aug(N_st+(5+kk)*N_st+3,ind);
%         pos_sen(kk,jj) = sqrt(x_sen(kk,jj)^2+y_sen(kk,jj)^2+z_sen(kk,jj)^2);
%         Hsen_px(kk,jj) = H_mat(1,6+kk,jj);
%         Hsen_py(kk,jj) = H_mat(2,6+kk,jj);
%         Hsen_pz(kk,jj) = H_mat(3,6+kk,jj);
%         Hsen_vx(kk,jj) = H_mat(4,6+kk,jj);
%         Hsen_vy(kk,jj) = H_mat(5,6+kk,jj);
%         Hsen_vz(kk,jj) = H_mat(6,6+kk,jj);
%         Hsen_p(kk,jj) = sqrt(Hsen_px(kk,jj)^2 + Hsen_py(kk,jj)^2 + Hsen_pz(kk,jj)^2);
%         Hsen_v(kk,jj) = sqrt(Hsen_vx(kk,jj)^2 + Hsen_vy(kk,jj)^2 + Hsen_vz(kk,jj)^2);
%     end
% end

for kk = 1:N_st     % Af_true for OFF
    for jj = 1:numel(time_pred_utc) %N_orb
        Hsen_px(kk,jj) = H_mat(1,kk,jj);
        Hsen_py(kk,jj) = H_mat(2,kk,jj);
        Hsen_pz(kk,jj) = H_mat(3,kk,jj);
        Hsen_vx(kk,jj) = H_mat(4,kk,jj);
        Hsen_vy(kk,jj) = H_mat(5,kk,jj);
        Hsen_vz(kk,jj) = H_mat(6,kk,jj);
        Hsen_p(kk,jj) = sqrt(Hsen_px(kk,jj)^2 + Hsen_py(kk,jj)^2 + Hsen_pz(kk,jj)^2);
        Hsen_v(kk,jj) = sqrt(Hsen_vx(kk,jj)^2 + Hsen_vy(kk,jj)^2 + Hsen_vz(kk,jj)^2);
    end
end

some_mat = H_mat;
X_state = [X_init;all_coeff'];
for jj = 1:numel(time_pred_utc)
    H_vec = [some_mat(1,:,jj)/sigma_pos;some_mat(2,:,jj)/sigma_pos;some_mat(3,:,jj)/sigma_pos;some_mat(4,:,jj)/sigma_vel;...
        some_mat(5,:,jj)/sigma_vel;some_mat(6,:,jj)/sigma_vel]';
    H_vec  = H_vec.*repmat(X_state,1,6);
    Hnorm_sen(:,jj) = vecnorm(H_vec,Inf,2);
    Hnorm_mat = repmat(vecnorm(H_vec,2,2),1,6);
    H_vec = H_vec./Hnorm_mat;
    for kk = 1:N_st
        Hdot(kk,:,jj) = H_vec(kk,:)*H_vec';
    end
end
Hd_mean = nanmean(abs(Hdot),3);

%% STM processing
% stm_sum = zeros(N_st,N_st);
% stm_old = eye(N_st);
% for jj = 1:numel(time_pred_utc)
%     stm_vec = X_aug(N_st+1:N_st+N_st^2,jj);
%     stm = reshape(stm_vec,N_st,N_st);
%     stm_sum = stm_sum + stm;
%     stm_mat = stm;
%     s_stm(:,jj) = svd(stm_mat);
%     s_stm_thresh(jj) = N_st*max(s_stm(:,jj))*2.220446049250313e-16;    
% %     stm_old = stm;
%     rank_stm(jj) = rank(stm_mat);
%     stm_state_vec = stm_mat(1:6,:)';
%     stmnorm(:,jj) = vecnorm(stm_state_vec,2,2);
%     stmnorm_mat = repmat(vecnorm(stm_state_vec,2,2),1,6);
%     stm_state_vec = stm_state_vec./stmnorm_mat;
%     for kk = 1:N_st
%         stmdotprod(kk,:,jj)= stm_state_vec(kk,:)*stm_state_vec';
%     end
%     stm_state_vec_mat(:,:,jj) = stm_state_vec;
% end
% stmd_mean = nanmean(abs(stmdotprod),3);
%%
% some_mat = H_mat;
% H_cumul = [];
% for jj = 1:numel(time_pred_utc)
%     H_cumul(6*jj-5:6*jj,:) = [some_mat(1,:,jj)/sigma_pos;some_mat(2,:,jj)/sigma_pos;some_mat(3,:,jj)/sigma_pos;...
%                               some_mat(4,:,jj)/sigma_vel;some_mat(5,:,jj)/sigma_vel;some_mat(6,:,jj)/sigma_vel];
%     H_vec = H_cumul';
%     Hnorm = repmat(vecnorm(H_vec,2,2),1,numel(H_vec(1,:)));
%     H_vec = H_vec./Hnorm;
% for kk = 1:N_st
%     Hdot(kk,:,jj) = H_vec(kk,:)*H_vec';
% end
% end

%%
% Hnorm_sen = stmnorm;
xvec = time_pred_utc/T_orb;
figure(2)
semilogy(xvec,Hnorm_sen(7,:),'b','LineWidth', 1)
hold on
semilogy(xvec,Hnorm_sen(8,:),'--b','LineWidth', 1)
semilogy(xvec,Hnorm_sen(9,:),'-.b','LineWidth', 1)
semilogy(xvec,Hnorm_sen(10,:),'k','LineWidth', 1)
semilogy(xvec,Hnorm_sen(11,:),'--k','LineWidth', 1)
semilogy(xvec,Hnorm_sen(12,:),'-.k','LineWidth', 1)
semilogy(xvec,Hnorm_sen(13,:),'r','LineWidth', 1)
semilogy(xvec,Hnorm_sen(14,:),'--r','LineWidth', 1)
% semilogy(xvec,Hnorm_sen(15,:),'-.r','LineWidth', 1)
% semilogy(xvec,Hnorm_sen(16,:),'c','LineWidth', 1)
% semilogy(xvec,Hnorm_sen(17,:),'m','LineWidth', 1)
semilogy(xvec, ones(1,numel(time_pred_utc)),'g','LineWidth',1)
title('Measurement sensitivity')
xlabel('No. of orbits')
ylabel('Sensitivity norm')
% leg = legend('$\overline{\mathcal{A}}_0$','$\overline{\mathcal{A}}_1$','$\overline{\mathcal{B}}_1$','$\overline{\mathcal{A}}_{2}$',...
%     '$\overline{\mathcal{A}}_{3}$','$\overline{\mathcal{B}}_{3}$','$\overline{\mathcal{A}}_{4}$','$\overline{\mathcal{A}}_{5}$','$\overline{\mathcal{B}}_{5}$'...
%     ,'$\overline{\mathcal{A}}_{6}$','$\overline{\mathcal{A}}_{7}$','$\overline{\mathcal{B}}_{7}$');
leg = legend('$\overline{\mathcal{A}}_0$','$\overline{\mathcal{A}}_4$','$\overline{\mathcal{A}}_8$','$\overline{\mathcal{A}}_{12}$',...
    '$\overline{\mathcal{A}}_{16}$','$\overline{\mathcal{A}}_{20}$','$\overline{\mathcal{A}}_{24}$','$\overline{\mathcal{A}}_{28}$');
% leg = legend('$\overline{{A}}_0$','$\overline{{A}}_1$','$\overline{{A}}_2$','$\overline{{A}}_3$',...
%     '$\overline{{A}}_4$','$\overline{{A}}_5$','$\overline{{A}}_6$','$\overline{{A}}_7$','$\overline{{A}}_8$',...
%     '$\overline{{A}}_9$','$\overline{{A}}_{10}$');
set(leg,'interpreter','latex')
set(gca,'FontSize',18)
xlim([0 17])
grid on

% GPM
% xvec = time_pred_utc/T_orb;
% figure(2)
% semilogy(xvec,Hnorm_sen(7,:),'b','LineWidth', 1)
% hold on
% semilogy(xvec,Hnorm_sen(8,:),'k','LineWidth', 1)
% semilogy(xvec,Hnorm_sen(9,:),'r','LineWidth', 1)
% semilogy(xvec,Hnorm_sen(10,:),'m','LineWidth', 1)
% semilogy(xvec,Hnorm_sen(11,:),'--b','LineWidth', 1)
% semilogy(xvec,Hnorm_sen(13,:),'--k','LineWidth', 1)
% semilogy(xvec,Hnorm_sen(29,:),'--r','LineWidth', 1)
% semilogy(xvec, ones(1,numel(time_pred_utc)),'g','LineWidth',1)
% title('Measurement sensitivity')
% xlabel('No. of orbits')
% ylabel('Sensitivity norm')
% leg = legend('$\overline{\mathcal{A}}_0$','$\overline{\mathcal{A}}_1$','$\overline{\mathcal{A}}_2$','$\overline{\mathcal{A}}_{3}$',...
%     '$\overline{\mathcal{A}}_{4}$','$\overline{\mathcal{A}}_{6}$','$\overline{\mathcal{B}}_{1}$');
% set(leg,'interpreter','latex')
% set(gca,'FontSize',18)
% xlim([0 17])
% grid on
%% Cube
% Hd_mean = stmd_mean;
% thresh = 0.99;
% [row_th,col_th] = find(Hd_mean(7:end,7:end)>thresh);
% figure()
% imagesc(Hd_mean(7:end,7:end))
% c = colorbar;
% c.Label.String = 'Dot product';
% title('Normalized dot product')
% xticks([1 2 3 4 5 6 7 8 9 10 11])
% yticks([1 2 3 4 5 6 7 8 9 10 11])
% set(gca,'TickLabelInterpreter', 'latex');
% set(gca,'XTickLabel', {'$\overline{A}_0$','$\overline{A}_1$','$\overline{A}_2$','$\overline{A}_3$',...
%     '$\overline{A}_4$','$\overline{A}_5$','$\overline{A}_6$','$\overline{A}_7$','$\overline{A}_8$',...
%     '$\overline{A}_9$','$\overline{A}_{10}$'});
% set(gca,'YTickLabel', {'$\overline{A}_0$','$\overline{A}_1$','$\overline{A}_2$','$\overline{A}_3$',...
%     '$\overline{A}_4$','$\overline{A}_5$','$\overline{A}_6$','$\overline{A}_7$','$\overline{A}_8$',...
%     '$\overline{A}_9$','$\overline{A}_{10}$'});
% set(gca,'FontSize',18)
% hold on
% plot(row_th,col_th,'r+', 'MarkerSize', 10);

%% GPM BODF
coeff_ind = 7:N_st;
ind_snr = Hnorm_sen(7:end,8641)>1;
coeff_ind = coeff_ind(ind_snr);

str_x = Xf_name_all(coeff_ind-6);
str_x = [{'$\overline{A}_{0}$'},{'$\overline{A}_{1}$'},{'$\overline{A}_{2}$'},{'$\overline{A}_{3}$'},{'$\overline{A}_{4}$'},...
             {'$\overline{A}_{5}$'}, {'$\overline{A}_{6}$'}, {'$\overline{A}_{7}$'}, {'$\overline{A}_{8}$'}, {'$\overline{A}_{9}$'}, {'$\overline{A}_{10}$'}];
% str_x = [{'$\overline{A}_{0,0}$'},{'$\overline{A}_{4,0}$'},{'$\overline{A}_{0,1}$'},{'$\overline{A}_{4,1}$'},{'$\overline{A}_{0,2}$'},{'$\overline{A}_{0,3}$'}];
% str_x = [{'\boldmath$A_{0,0}$'},{'\boldmath$A_{1,0}$'},{'\boldmath$A_{2,0}$'},{'\boldmath$A_{3,0}$'},{'\boldmath$A_{4,0}$'},...
%     {'\boldmath$A_{6,0}$'},{'\boldmath$A_{0,1}$'},{'\boldmath$A_{1,1}$'},{'\boldmath$A_{2,1}$'},{'\boldmath$A_{3,1}$'},{'\boldmath$A_{4,1}$'},...
%     {'\boldmath$A_{0,2}$'},{'\boldmath$A_{1,2}$'},{'\boldmath$A_{2,2}$'},{'\boldmath$A_{3,2}$'},{'\boldmath$A_{0,3}$'},{'\boldmath$A_{1,3}$'},...
%     {'\boldmath$A_{2,3}$'},{'\boldmath$A_{4,3}$'},{'\boldmath$A_{2,4}$'},{'\boldmath$A_{2,5}$'},{'\boldmath$B_{1,0}$'},{'\boldmath$B_{1,1}$'}];

%%%%%%%%
figure()
h = bar(Hnorm_sen(coeff_ind,end))
set(gca,'yscale','log')
set(gca,'FontSize',18)
ylabel('Sensitivity norm')
title('Measurement sensitivities')
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTickLabel', str_x);
% yticks([1e-4, 1e-2, 1e0, 1e2, 1e4,1e6])
xticks([1:numel(coeff_ind)])
xtickangle(45)
% ylim([1e-10 1e7])
grid on
ax = ancestor(h, 'axes');
Xrule = ax.XAxis;
% Xrule.FontSize = 12;