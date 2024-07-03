addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Density inversion methods')
clc
clearvars
load champ_jb08_g5
X_jb08 = X_true_aug;
rho_jb08 = rho_avg;
load champ_msis_g5
X_msis = X_true_aug;
rho_msis = rho_avg;
load champ_hasdm_g5
X_hasdm = X_true_aug;
rho_hasdm_calc = rho_avg;
load champ_g5
X_champ = X_true_aug;
rho_champ = rho_avg;

doy_g5 = doy;
err_msis = X_msis - X_champ;
err_jb08 = X_jb08 - X_champ;
err_nodrag = X_hasdm - X_champ(:,1:numel(X_hasdm(1,:))) ;
X_ref = X_hasdm;
for ii = 1:numel(time_prop_utc)
    rhat = X_ref(1:3,ii)/norm(X_ref(1:3,ii));
vhat = X_ref(4:6,ii)/norm(X_ref(4:6,ii));
    Xhat = rhat;
Zhat = cross(rhat,vhat);
Yhat = cross(Zhat,Xhat);
Yhat = Yhat/norm(Yhat);
ECI2RIC = [Xhat';Yhat';Zhat'];
err_msis_ric(:,ii) = ECI2RIC*err_msis(1:3,ii);
err_jb08_ric(:,ii) = ECI2RIC*err_jb08(1:3,ii);
err_nodrag_ric(:,ii) = ECI2RIC*err_nodrag(1:3,ii);

coe = rv2coe_E(X_champ(1:3,ii),X_champ(4:6,ii),mu_e);
a = coe(1);
n_champ = sqrt(mu_e/a^3); 
v_champ = norm(X_champ(4:6,ii));
M_champ(ii) = n_champ*time_prop_utc(ii);

coe = rv2coe_E(X_hasdm(1:3,ii),X_hasdm(4:6,ii),mu_e);
a = coe(1);
n_hasdm = sqrt(mu_e/a^3); 
v_hasdm = norm(X_hasdm(4:6,ii));
M_hasdm(ii) = n_hasdm*time_prop_utc(ii);

coe = rv2coe_E(X_msis(1:3,ii),X_msis(4:6,ii),mu_e);
a = coe(1);
n_msis = sqrt(mu_e/a^3); 
M_msis(ii) = n_msis*time_prop_utc(ii);

coe = rv2coe_E(X_jb08(1:3,ii),X_jb08(4:6,ii),mu_e);
a = coe(1);
n_jb08 = sqrt(mu_e/a^3); 
M_jb08(ii) = n_jb08*time_prop_utc(ii);

del_M = (M_champ(ii) - M_msis(ii));
del_s_msis(ii) = del_M*v_champ/n_champ;

del_M = (M_champ(ii) - M_jb08(ii));
del_s_jb08(ii) = del_M*v_champ/n_champ;


del_M = (M_champ(ii) - M_hasdm(ii));
del_s_hasdm(ii) = del_M*v_champ/n_champ;
end
%%
figure(1)
subplot(2,1,1)
plot(time_avg/86400,rho_msis,'LineWidth',2)
hold on
plot(time_avg/86400,rho_jb08,'LineWidth',2)
plot(time_avg/86400,rho_hasdm_calc,'LineWidth',2)
plot(time_avg/86400,rho_champ,'k','LineWidth',4)
grid on
legend('NRLMSISE-00','JB2008','HASDM','CHAMP')
title('Density (kg/m^3)');
set(gca,'FontSize',16)
subplot(2,1,2)
plot(time_prop_utc/86400,del_s_msis/1e3,'LineWidth',2)
hold on
plot(time_prop_utc/86400,del_s_jb08/1e3,'LineWidth',2)
plot(time_prop_utc/86400,del_s_hasdm/1e3,'LineWidth',2)
grid on
% legend('NRLMSISE-00','JB2008','HASDM')
title('In-track error w.r.t CHAMP density (km)');
xlabel('Days')
set(gca,'FontSize',16)


