%% Kaula's rule for OFF coefficients
clc
clearvars
load('Cd_sphere_dria_150F10','Af_total_mat','Af_mean','Af_std')

Af_max = max(abs(Af_total_mat),[],2);
Y1_mean = abs(Af_mean(2:2:end));
Y2_mean = abs(Af_mean(3:2:end)); 
Y1_std = abs(Af_mean(2:2:end))+abs(Af_std(2:2:end));
Y2_std = abs(Af_mean(3:2:end))+abs(Af_std(3:2:end));

Y1_upper = abs(Af_max(2:2:end));
Y2_upper = abs(Af_max(3:2:end));
N = numel(Af_total_mat(1,:));
X1 = 1:5;
X2 = 1:5;
model_fun = @(a,x)(exp(a(1)*x+a(2)));  %

a0 = [-1,-1];
a1_mean = nlinfit(X1,Y1_mean',model_fun,a0);
a2_mean = nlinfit(X2,Y2_mean',model_fun,a0);

a1_std = nlinfit(X1,Y1_std',model_fun,a0);
a2_std = nlinfit(X2,Y2_std',model_fun,a0);

a1_upper = nlinfit(X1,Y1_upper',model_fun,a0);
a2_upper = nlinfit(X2,Y2_upper',model_fun,a0);

%%
figure(1)
semilogy(model_fun(a1_mean,X1),'r','LineWidth',1)
hold on
semilogy(ones(N),abs(Af_total_mat(2,:)),'.b','LineWidth',1)
semilogy(2*ones(N),abs(Af_total_mat(4,:)),'.b','LineWidth',1)
semilogy(3*ones(N),abs(Af_total_mat(6,:)),'.b','LineWidth',1)
semilogy(4*ones(N),abs(Af_total_mat(8,:)),'.b','LineWidth',1)
semilogy(5*ones(N),abs(Af_total_mat(10,:)),'.b','LineWidth',1)
leg = legend('$e^{-(2.9k+0.3)}$','Distribution');
set(leg,'interpreter','latex')
xlabel('k')
ylabel('Magnitude')
title('Amplitude of even-order Fourier coefficients')
set(gca,'FontSize',18)

figure(2)
semilogy(model_fun(a2_mean,X2),'r','LineWidth',1)
hold on
semilogy(ones(N),abs(Af_total_mat(3,:)),'.b','LineWidth',1)
semilogy(2*ones(N),abs(Af_total_mat(5,:)),'.b','LineWidth',1)
semilogy(3*ones(N),abs(Af_total_mat(7,:)),'.b','LineWidth',1)
semilogy(4*ones(N),abs(Af_total_mat(9,:)),'.b','LineWidth',1)
semilogy(5*ones(N),abs(Af_total_mat(11,:)),'.b','LineWidth',1)
leg = legend('$e^{-(3.9k+0.6)}$','Distribution');
set(leg,'interpreter','latex')
xlabel('k')
ylabel('Magnitude')
title('Amplitude of odd-order Fourier coefficients')
set(gca,'FontSize',18)