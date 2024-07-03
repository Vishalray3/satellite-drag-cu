%% ISS plots
ss = mod(time_prop_utc, 60);
min = mod(floor(time_prop_utc/60),60);
hr = mod(floor(time_prop_utc/3600),24);
dy = day_mon+floor(time_prop_utc/86400); 
mn = mon*ones(1,numel(dy));
yr = yyyy*ones(1,numel(dy));

t_iss = datetime(yr,mn,dy,hr,min,ss);

%%
figure()
load iss_trajectory_dec14-16_2006
plot(t_iss,rho, 'k','LineWidth', 1)
hold on
load iss_trajectory_dec14-16_2006_jb08
plot(t_iss,rho, 'b','LineWidth', 1)
load iss_trajectory_dec14-16_2006_msis
plot(t_iss,rho, 'r','LineWidth', 1)
ylabel('Density ($kg/m^3$)','interpreter','latex')
legend('HASDM','JB2008','NRLMSISE-00')
title('Atmospheric density in ISS orbit')
set(gca,'FontSize',16)
grid on