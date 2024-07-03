%% Torque Fourier coefficient estimation
rng('default')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Density inversion methods')
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Density inversion methods/Simultaneous estimation/Torque estimation')
clc
clearvars
Hp_ind = 305;
Ha_ind = 310;
Constants_torque
load('torque_fourier_case1')
load('300km_torque_case1', 'reci', 'veci','rho','theta','time_prop_utc','omega_e','rho_msis')
load('300km_torquefourier_noise')
rho_est = rho_msis;
theta = theta*180/pi;
Order_b = 30;
Pf_std = 0.2*abs(Pf_total);
Qf_std = 0.2*abs(Qf_total);
R_aug = diag([1e-10,1e-10,1e-10]);
% Pf_noise = normrnd(0,Pf_std);
% Qf_noise = normrnd(0,Qf_std);
Pf_total_err = Pf_total + Pf_noise;
Qf_total_err = Qf_total + Qf_noise;
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
str_est1 = [{'Q11'},{'Q13'},{'Q15'}];
str_est2 = [{'P21'},{'P23'},{'P25'}];
str_est3 = [{'P31'},{'P32'},{'P33'},{'P35'},{'P37'},{'Q31'}];

vec_est1 = find(ismember(Xf_name1,str_est1));    N_est1 = numel(vec_est1);
vec_est2 = find(ismember(Xf_name2,str_est2));    N_est2 = numel(vec_est2);
vec_est3 = find(ismember(Xf_name3,str_est3));    N_est3 = numel(vec_est3);

P_prior = diag([Xf_std1(vec_est1).^2;Xf_std2(vec_est2).^2;Xf_std3(vec_est3).^2]);
Hmat = [];
Ymat = [];
x_prior = zeros(N_est1+N_est2+N_est3,1);
x_est = x_prior;
iter = 10;
N_plates = numel(An_total(1,:));
ovec_cd = [0:Order_b+1]';
ovec = [0:Order_b];

err_init = [X_f1t(vec_est1); X_f2t(vec_est2); X_f3t(vec_est3)] - [X_f1(vec_est1); X_f2(vec_est2); X_f3(vec_est3)];


for n = 1:iter
    x_prior = x_prior - x_est;
    P_low = chol(P_prior, 'lower');
    Lambda = inv(P_low')*inv(P_low);
    N =   Lambda*x_prior;
    
    X_est1 = Xf_est(1:N_f1); X_est2 = Xf_est(N_f1+1:N_f1+N_f2); X_est3 = Xf_est(N_f1+N_f2+1:end);
    X_est1(vec_est1) = X_est1(vec_est1) + x_est(1:N_est1);
    X_est2(vec_est2) = X_est2(vec_est2) + x_est(N_est1+1:N_est1+N_est2);
    X_est3(vec_est3) = X_est3(vec_est3) + x_est(N_est1+N_est2+1:N_est1+N_est2+N_est3);
    Xf_est = [X_est1;X_est2;X_est3];
    for ii = 1:numel(time_prop_utc)
        V_vec = [cosd(theta(ii)); sind(theta(ii)); 0];
        V_mat = repmat(V_vec, 1, N_plates);
        Cd_plates = sum(An_total.*repmat(cosd(ovec_cd*theta(ii)),1,N_plates),1) + sum(Bn_total.*repmat(sind(ovec_cd*theta(ii)),1,N_plates),1);
        
        v_r = veci(:,ii) - cross([0;0;omega_e],reci(:,ii));
        Vi = norm(v_r);
        
        torque_actual(:,ii) = -0.5*rho(ii)*Vi^2/(mass)*sum(cross(Cp, V_mat, 1).*repmat(Cd_plates,3,1) ,2);
        
        Xf_arg = [cosd(ovec*theta(ii))';sind(ovec*theta(ii))'];
        Xf_arg1 = Xf_arg(Xf_ind1);   Xf_arg2 = Xf_arg(Xf_ind2);    Xf_arg3 = Xf_arg(Xf_ind3);
        
        X_est1 = Xf_est(1:N_f1); X_est2 = Xf_est(N_f1+1:N_f1+N_f2); X_est3 = Xf_est(N_f1+N_f2+1:end);
        Tx = sum(X_est1.*Xf_arg1); Ty = sum(X_est2.*Xf_arg2); Tz = sum(X_est3.*Xf_arg3);
        
        torque_predict(:,ii) =  -0.5*rho_est(ii)*Vi^2/(mass)*[Tx;Ty;Tz];
        
        Y_tor(:,ii) = torque_actual(:,ii) - torque_predict(:,ii);
        residual = Y_tor(:,ii);
        
        H = -0.5*rho_est(ii)*Vi^2/(mass)*[Xf_arg1(vec_est1)', zeros(1,N_est2), zeros(1,N_est3); zeros(1,N_est1), Xf_arg2(vec_est2)', zeros(1,N_est3);...
            zeros(1,N_est1),zeros(1,N_est2),Xf_arg3(vec_est3)'];
        
        Lambda = Lambda + H'/R_aug*H;
        Lambda_low = chol(Lambda,'lower');
        Lambda_inv = inv(Lambda_low')*inv(Lambda_low);
        P_up = Lambda_inv;
        sigma_x(:,ii) = sqrt(diag(P_up));
        
        
        N = N + H'/R_aug*residual;
        
        Xest_mat(:,ii) = Xf_est;
        residual_mat(n).residual(:,ii) = residual;
    end
    x_est = Lambda_inv*N;
    
    x_est
    if norm(x_est) < 5e-3
        break
    end
    
end

err_final = [X_f1t(vec_est1); X_f2t(vec_est2); X_f3t(vec_est3)] - [X_est1(vec_est1); X_est2(vec_est2); X_est3(vec_est3)];

%% Plots
figure(1)
subplot(3,1,1)
plot(time_prop_utc/3600,residual_mat(1).residual(1,:), 'LineWidth',2)
hold on
plot(time_prop_utc/3600,residual_mat(n).residual(1,:), 'LineWidth',2)
ylabel('X (Nm)')
title('Torque estimation errors')
legend('Initial','Estimation')
set(gca,'FontSize',16)
subplot(3,1,2)
plot(time_prop_utc/3600,residual_mat(1).residual(2,:), 'LineWidth',2)
hold on
plot(time_prop_utc/3600,residual_mat(n).residual(2,:), 'LineWidth',2)
ylabel('Y (Nm)')
set(gca,'FontSize',16)
subplot(3,1,3)
plot(time_prop_utc/3600,residual_mat(1).residual(3,:), 'LineWidth',2)
hold on
plot(time_prop_utc/3600,residual_mat(n).residual(3,:), 'LineWidth',2)
ylabel('Z (Nm)')
xlabel('Time (hours)')
set(gca,'FontSize',16)

figure(2)
plot(time_prop_utc/3600,rho, 'LineWidth',2)
hold on
plot(time_prop_utc/3600,rho_msis, 'LineWidth',2)
ylabel('Density ($kg/m^3$)','Interpreter','latex')
title('Density')
legend('True (HASDM)','Filter baseline (NRLMSISE-00)')
grid on
set(gca,'FontSize',16)
%%
% vec_ord = [0:Order_b];
% figure(3)
% subplot(2,1,1)
% plot(vec_ord, Pf_total(1,:), 'LineWidth',2)
% hold on
% plot(vec_ord, Pf_total(2,:), 'LineWidth',2)
% plot(vec_ord, Pf_total(3,:), 'LineWidth',2)
% ylabel('Cosine coefficients')
% set(gca,'FontSize',16)
% title('Torque Fourier coefficients')
% legend('x','y','z')
% grid on
% subplot(2,1,2)
% plot(vec_ord, Qf_total(1,:), 'LineWidth',2)
% hold on
% plot(vec_ord, Qf_total(2,:), 'LineWidth',2)
% plot(vec_ord, Qf_total(3,:), 'LineWidth',2)
% ylabel('Sine coefficients')
% xlabel('Order')
% set(gca,'FontSize',16)
% grid on