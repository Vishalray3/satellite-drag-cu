%% Fitting an autocovariance function
clc
clearvars
close all
% load('2_batch_msis_cdcr','rho_est','rho_true','T_orb','time_prop_utc')
load('spire83_propagation_2018_11_7','rho_hasdm','rho_msis')
% load('GPM_data2_Eric_hasdm','hasden')
% generate errors for the true density
% t_ang = time_prop_utc/T_orb*2*pi; 
% rho_err = 1e-12*sin(t_ang);
rho_err = rho_msis-rho_hasdm; % 0.1*rho_true;
% rho_err = rho_err;

% compute the autocovariance
maxlag = 2000;
[autocov, lags] = xcov(rho_err,maxlag,'normalized');
lags = lags*10;
[autocov2, ~] = xcov(rho_err,maxlag);

%% fit the autocovariance function
fo = fitoptions('Method','NonlinearLeastSquares','Normalize','off',...
    'Lower',[0,0,0,0,0],...
    'Upper',[Inf,1,Inf, Inf,Inf],...
    'StartPoint',[1,0,1e-3, 0.003, 0.004]);
ft = fittype('a^2*exp(-b*c*abs(x))*(cos(c*sqrt(1-b^2)*abs(x)) + b*c/(c*sqrt(1-b^2))*sin(c*sqrt(1-b^2)*abs(x))) + d^2*exp(-abs(x)*e)*sqrt(pi)/2/e','options',fo);
% fo = fitoptions('Method','NonlinearLeastSquares',...
%     'Lower',[0,0,0],...
%     'Upper',[Inf,1,Inf],...
%     'StartPoint',[1,0,0.1]);
% ft = fittype('a^2*exp(-b*c*abs(x))*(cos(c*sqrt(1-b^2)*abs(x)) + b*c/(c*sqrt(1-b^2))*sin(c*sqrt(1-b^2)*abs(x)))','options',fo);
% ft = fittype('a^2*exp(-c*d*abs(x))*cos(d*sqrt(1-c^2)*abs(x))','options',fo);
% fo = fitoptions('Method','NonlinearLeastSquares',...
%     'Lower',[-Inf,-Inf],...
%     'Upper',[Inf,Inf],...
%     'StartPoint',[1,1]);
% ft = fittype('a^2*cosd(c*x)','options',fo);
% x = a^2*exp(-b*c*abs(x))*(cos(c*sqrt(1-b^2)*abs(x)) + b*c/(c*sqrt(1-b^2))*sin(c*sqrt(1-b^2)*abs(x)));
ind = find(lags == 0);
nor_fac = autocov2(ind);
[curve2,gof2] = fit(lags',autocov',ft);
sigma = curve2.a*sqrt(nor_fac)
zeta = curve2.b
omega_n = curve2.c
c2 = curve2.d*sqrt(nor_fac)
tau_inv_x3 = curve2.e
chi = 0; %curve2.b;
% alpha = atan(zeta/sqrt(1-zeta^2));
% c1 = sqrt(2*sigma^2*omega_n*sin(alpha-chi)/cos(chi))
% c = sqrt(2*sigma^2*omega_n^3*sin(alpha+chi)/cos(chi));
% c2 = c - 2*c1*zeta*omega_n
c= sqrt(4*omega_n^3*zeta*sigma^2)
%%
curve2
plot(lags,autocov)
hold on
plot(lags,curve2(lags))

% figure(2)
% plot(lags,autocov2)
% hold on
% plot(lags,nor_fac*curve2(lags))
% sig = 1e-10;
% zeta = 0;
% omega_n = 0.02;
% beta = omega_n*sqrt(1-zeta^2);
% plot(sig^2*exp(-zeta*omega_n*abs(lags)).*(cos(beta*abs(lags)) + zeta*omega_n/beta*sin(beta*abs(lags))))

%%
% 
% ind = find(lags == 0);
% nor_fac = autocov2(ind);
% sigma = sqrt(nor_fac);
% zeta = 0.1;
% omega_n = 2*pi/T_orb*2;
% cov = auto_cov(sqrt(nor_fac),0.1, 2*pi/T_orb*2, lags);
% c= sqrt(4*omega_n^3*zeta*sigma^2)
% plot(autocov2)
% hold on
% plot(cov)
% function [cov] = auto_cov(sigma, zeta, omega_n,lags)
% cov = sigma^2*exp(-zeta*omega_n*abs(lags)).*(cos(omega_n*sqrt(1-zeta^2)*abs(lags)) + zeta*omega_n/(omega_n*sqrt(1-zeta^2))*sin(omega_n*sqrt(1-zeta^2)*abs(lags)));
% end
