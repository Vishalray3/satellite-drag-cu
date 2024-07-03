%% Fitting an autocovariance function
clc
clearvars
load Truth_msis_circ_gpm_tseries

% generate errors for the true density
t_ang = time_prop_utc/T_orb*2*pi; 

% compute the autocovariance
maxlag = 100;
[autocov, lags] = xcov(Bf_tseries(2,:),maxlag,'normalized');
lags = lags*10;
[autocov2, ~] = xcov(Bf_tseries(2,:),maxlag);

%% fit the autocovariance function
fo = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0, 0,0],...
    'Upper',[Inf,1,Inf],...
    'StartPoint',[1,0,0.1]);
ft = fittype('a^2*exp(-b*c*abs(x))*(cos(c*sqrt(1-b^2)*abs(x)) + b*c/(c*sqrt(1-b^2))*sin(c*sqrt(1-b^2)*abs(x)))','options',fo);
% b = zeta, c = omega
[curve2,gof2] = fit(lags',autocov',ft);
%%
ind = find(lags == 0);
nor_fac = autocov2(ind);
% sigma = sqrt(nor_fac)
plot(lags,autocov)
hold on
plot(lags,curve2(lags))

figure(2)
plot(lags,autocov2)
hold on
plot(lags,nor_fac*curve2(lags))
curve2
nor_fac
% sig = 1e-10;
% zeta = 0;
% omega_n = 0.02;
% beta = omega_n*sqrt(1-zeta^2);
% plot(sig^2*exp(-zeta*omega_n*abs(lags)).*(cos(beta*abs(lags)) + zeta*omega_n/beta*sin(beta*abs(lags))))