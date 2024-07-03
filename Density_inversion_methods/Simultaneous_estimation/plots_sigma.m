%% plots 
close all
vec = 1:numel(sigma_x(1,:));
figure()
% plot(time_prop_utc(vec)/3600,sigma_x(7,:),'LineWidth', 2)
% hold on
plot(time_prop_utc(vec)/3600,sigma_x(8,:),'LineWidth', 2)
hold on
plot(time_prop_utc(vec)/3600,sigma_x(9,:),'LineWidth', 2)
plot(time_prop_utc(vec)/3600,sigma_x(10,:),'LineWidth', 2)
plot(time_prop_utc(vec)/3600,sigma_x(11,:),'LineWidth', 2)
plot(time_prop_utc(vec)/3600,sigma_x(18,:),'LineWidth', 2)
plot(time_prop_utc(vec)/3600,sigma_x(19,:),'LineWidth', 2)
plot(time_prop_utc(vec)/3600,sigma_x(20,:),'LineWidth', 2)
plot(time_prop_utc(vec)/3600,sigma_x(21,:),'LineWidth', 2)
xlabel('Time (hours)')
ylabel('1-sigma')
legend('A1', 'A2', 'A3','A4','B1','B2','B3','B4')
title('Uncertainty of Fourier coefficients')
set(gca,'FontSize',18)
grid on
% xlim([0 23])