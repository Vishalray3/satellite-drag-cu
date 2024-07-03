%% Starlink simulation in-track error calculation
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Density inversion methods')
clc
clearvars
load starlink_250km_jb08_openbook
X_jb08 = X_true_aug;
load starlink_250km_msis_openbook
X_msis = X_true_aug;
load starlink_210km_nodrag
X_nodrag = X_true_aug;
load starlink_210km_hasdm
X_hasdm = X_true_aug;

doy_g5 = doy;
% err_msis = X_msis - X_hasdm;
% err_jb08 = X_jb08 - X_hasdm;
err_nodrag = X_nodrag(:,1:numel(X_hasdm(1,:))) - X_hasdm;
X_ref = X_hasdm;
for ii = 1:numel(time_prop_utc)
    rhat = X_ref(1:3,ii)/norm(X_ref(1:3,ii));
vhat = X_ref(4:6,ii)/norm(X_ref(4:6,ii));
    Xhat = rhat;
Zhat = cross(rhat,vhat);
Yhat = cross(Zhat,Xhat);
Yhat = Yhat/norm(Yhat);
ECI2RIC = [Xhat';Yhat';Zhat'];
% err_msis_ric(:,ii) = ECI2RIC*err_msis(1:3,ii);
% err_jb08_ric(:,ii) = ECI2RIC*err_jb08(1:3,ii);
err_nodrag_ric(:,ii) = ECI2RIC*err_nodrag(1:3,ii);

coe = rv2coe_E(X_hasdm(1:3,ii),X_hasdm(4:6,ii),mu_e);
a = coe(1);
n_hasdm = sqrt(mu_e/a^3); 
v_hasdm = norm(X_hasdm(4:6,ii));
M_hasdm(ii) = n_hasdm*time_prop_utc(ii);

coe = rv2coe_E(X_nodrag(1:3,ii),X_nodrag(4:6,ii),mu_e);
a = coe(1);

n_nodrag = sqrt(mu_e/a^3); 
v_nodrag = norm(X_nodrag(4:6,ii));
M_nodrag(ii) = n_nodrag*time_prop_utc(ii);

del_M(ii) = (M_hasdm(ii) - M_nodrag(ii));
del_s(ii) = del_M(ii)*v_nodrag/n_nodrag;
end

%%
load starlink_250km_jb08_openbook
X_jb08 = X_true_aug;
load starlink_250km_msis_openbook
X_msis = X_true_aug;
load starlink_210km_nodrag
X_nodrag = X_true_aug;
load starlink_210km_hasdm_g1
X_hasdm = X_true_aug;
doy_g1 = doy;
% err_msis = X_msis - X_hasdm;
% err_jb08 = X_jb08 - X_hasdm;
err_nodrag = X_nodrag(:,1:numel(X_hasdm(1,:))) - X_hasdm;
X_ref = X_hasdm;
for ii = 1:numel(time_prop_utc)
    rhat = X_ref(1:3,ii)/norm(X_ref(1:3,ii));
vhat = X_ref(4:6,ii)/norm(X_ref(4:6,ii));
    Xhat = rhat;
Zhat = cross(rhat,vhat);
Yhat = cross(Zhat,Xhat);
Yhat = Yhat/norm(Yhat);
ECI2RIC = [Xhat';Yhat';Zhat'];
% err_msis_ric_g1(:,ii) = ECI2RIC*err_msis(1:3,ii);
% err_jb08_ric_g1(:,ii) = ECI2RIC*err_jb08(1:3,ii);
err_nodrag_ric_g1(:,ii) = ECI2RIC*err_nodrag(1:3,ii);

coe = rv2coe_E(X_hasdm(1:3,ii),X_hasdm(4:6,ii),mu_e);
a = coe(1);
n_hasdm = sqrt(mu_e/a^3); 
v_hasdm = norm(X_hasdm(4:6,ii));
M_hasdm(ii) = n_hasdm*time_prop_utc(ii);

coe = rv2coe_E(X_nodrag(1:3,ii),X_nodrag(4:6,ii),mu_e);
a = coe(1);

n_nodrag = sqrt(mu_e/a^3); 
v_nodrag = norm(X_nodrag(4:6,ii));
M_nodrag(ii) = n_nodrag*time_prop_utc(ii);

del_M(ii) = (M_hasdm(ii) - M_nodrag(ii));
del_s_g1(ii) = del_M(ii)*v_nodrag/n_nodrag;
end

%%
time_vec = [0:40]/10;
load NRLMSISE_2003
load('starlink_210km_hasdm', 'doy')
doy_g5 = doy;
figure(2)
subplot(3,1,1)
load starlink_210km_hasdm
plot(time_avg/86400,r_alt_min/1e3,'k','LineWidth',2)
hold on
% load starlink_250km_msis_openbook
% plot(time_avg/86400,r_alt_avg/1e3,'k','LineWidth',2)
% load starlink_250km_jb08_openbook
% plot(time_avg/86400,r_alt_avg/1e3,'r','LineWidth',2)
% load starlink_250km_hasdm_openbook
% plot(time_avg/86400,r_alt_avg/1e3,'-g','LineWidth',2)
load starlink_210km_hasdm_g1
plot(time_avg/86400,r_alt_min/1e3,'--k','LineWidth',2)
% load starlink_250km_msis_openbook_g1
% plot(time_avg/86400,r_alt_avg/1e3,'--k','LineWidth',2)
% load starlink_250km_jb08_openbook_g1
% plot(time_avg/86400,r_alt_avg/1e3,'--r','LineWidth',2)
grid on
title('Perigee altitude decay (km)')
legend('G5 storm','G1 storm')
set(gca,'FontSize',14)
yticks([170, 190, 210])
ylim([170 210])
subplot(3,1,2)
plot(time_prop_utc/86400,del_s/1e3,'k','LineWidth',2)
hold on
plot(time_prop_utc/86400,del_s_g1/1e3,'--k','LineWidth',2)
% plot(time_prop_utc/86400,err_msis_ric_g1(2,:)/1e3,'--k','LineWidth',2)
% plot(time_prop_utc/86400,err_jb08_ric_g1(2,:)/1e3,'--r','LineWidth',2)
grid on

title('In-track position difference (km)')
% legend('G5 storm','G1 storm')
set(gca,'FontSize',14)
subplot(3,1,3)
yyaxis left
plot(time_vec,Ap_total((doy_g5-1)*8:8*(doy_g5+4)),'k-','LineWidth',2)
hold on
plot(time_vec,Ap_total((doy_g1-1)*8:8*(doy_g1+4)),'k--','LineWidth',2)
% legend('G5 storm','G1 storm')
set(gca,'YColor', 'k')
ylabel('Daily Ap')
yyaxis right
plot([0:4],F10_total(doy_g5:doy_g5+4),'-','LineWidth',2)
hold on
plot([0:4],F10_total(doy_g1:doy_g1+4),'--','LineWidth',2)
xlabel('Days')
ylabel('Daily F10.7')
title('Space weather indices')
grid on
set(gca,'FontSize',14)


%%
clc
clearvars
vec = 1:8641;
load starlink_250km_hasdm_openbook
X_hasdm = X_true_aug(:,vec);
load starlink_250km_jb08_openbook
X_jb08 = X_true_aug(:,vec);
load starlink_250km_msis_openbook
X_msis = X_true_aug(:,vec);
doy_g5 = doy;
err_msis = X_msis - X_hasdm;
err_jb08 = X_jb08 - X_hasdm;
for ii = 1:numel(vec)


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

del_M = (M_hasdm(ii) - M_msis(ii));
del_s_msis_250(ii) = del_M*v_hasdm/n_hasdm;

del_M = (M_hasdm(ii) - M_jb08(ii));
del_s_jb08_250(ii) = del_M*v_hasdm/n_hasdm;
end

load starlink_250km_hasdm_openbook_g1
X_hasdm = X_true_aug(:,vec);
load starlink_250km_jb08_openbook_g1
X_jb08 = X_true_aug(:,vec);
load starlink_250km_msis_openbook_g1
X_msis = X_true_aug(:,vec);
doy_g5 = doy;
err_msis = X_msis - X_hasdm;
err_jb08 = X_jb08 - X_hasdm;
for ii = 1:numel(vec)


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

del_M = (M_hasdm(ii) - M_msis(ii));
del_s_msis_g1_250(ii) = del_M*v_hasdm/n_hasdm;

del_M = (M_hasdm(ii) - M_jb08(ii));
del_s_jb08_g1_250(ii) = del_M*v_hasdm/n_hasdm;
end

load starlink_350km_hasdm_openbook
X_hasdm = X_true_aug(:,vec);
load starlink_350km_jb08_openbook
X_jb08 = X_true_aug(:,vec);
load starlink_350km_msis_openbook
X_msis = X_true_aug(:,vec);
doy_g5 = doy;
err_msis = X_msis - X_hasdm;
err_jb08 = X_jb08 - X_hasdm;
for ii = 1:numel(vec)
    rhat = X_hasdm(1:3,ii)/norm(X_hasdm(1:3,ii));
vhat = X_hasdm(4:6,ii)/norm(X_hasdm(4:6,ii));
    Xhat = rhat;
Zhat = cross(rhat,vhat);
Yhat = cross(Zhat,Xhat);
Yhat = Yhat/norm(Yhat);
ECI2RIC = [Xhat';Yhat';Zhat'];
err_msis_ric_350(:,ii) = ECI2RIC*err_msis(1:3,ii);
err_jb08_ric_350(:,ii) = ECI2RIC*err_jb08(1:3,ii);

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

del_M = (M_hasdm(ii) - M_msis(ii));
del_s_msis(ii) = del_M*v_hasdm/n_hasdm;

del_M = (M_hasdm(ii) - M_jb08(ii));
del_s_jb08(ii) = del_M*v_hasdm/n_hasdm;
end

load starlink_350km_hasdm_openbook_g1
X_hasdm = X_true_aug(:,vec);
load starlink_350km_jb08_openbook_g1
X_jb08 = X_true_aug(:,vec);
load starlink_350km_msis_openbook_g1
X_msis = X_true_aug(:,vec);
doy_g5 = doy;
err_msis = X_msis - X_hasdm;
err_jb08 = X_jb08 - X_hasdm;
for ii = 1:numel(vec)
    rhat = X_hasdm(1:3,ii)/norm(X_hasdm(1:3,ii));
vhat = X_hasdm(4:6,ii)/norm(X_hasdm(4:6,ii));
    Xhat = rhat;
Zhat = cross(rhat,vhat);
Yhat = cross(Zhat,Xhat);
Yhat = Yhat/norm(Yhat);
ECI2RIC = [Xhat';Yhat';Zhat'];
err_msis_ric_350_g1(:,ii) = ECI2RIC*err_msis(1:3,ii);
err_jb08_ric_350_g1(:,ii) = ECI2RIC*err_jb08(1:3,ii);

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

del_M = (M_hasdm(ii) - M_msis(ii));
del_s_msis_g1(ii) = del_M*v_hasdm/n_hasdm;

del_M = (M_hasdm(ii) - M_jb08(ii));
del_s_jb08_g1(ii) = del_M*v_hasdm/n_hasdm;
end

load starlink_550km_hasdm_openbook
X_hasdm = X_true_aug(:,vec);
load starlink_550km_jb08_openbook
X_jb08 = X_true_aug(:,vec);
load starlink_550km_msis_openbook
X_msis = X_true_aug(:,vec);
doy_g5 = doy;
err_msis = X_msis - X_hasdm;
err_jb08 = X_jb08 - X_hasdm;
for ii = 1:numel(vec)
    rhat = X_hasdm(1:3,ii)/norm(X_hasdm(1:3,ii));
vhat = X_hasdm(4:6,ii)/norm(X_hasdm(4:6,ii));
    Xhat = rhat;
Zhat = cross(rhat,vhat);
Yhat = cross(Zhat,Xhat);
Yhat = Yhat/norm(Yhat);
ECI2RIC = [Xhat';Yhat';Zhat'];
err_msis_ric_550(:,ii) = ECI2RIC*err_msis(1:3,ii);
err_jb08_ric_550(:,ii) = ECI2RIC*err_jb08(1:3,ii);

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

del_M = (M_hasdm(ii) - M_msis(ii));
del_s_msis_550(ii) = del_M*v_hasdm/n_hasdm;

del_M = (M_hasdm(ii) - M_jb08(ii));
del_s_jb08_550(ii) = del_M*v_hasdm/n_hasdm;
end

load starlink_550km_hasdm_openbook_g1
X_hasdm = X_true_aug(:,vec);
load starlink_550km_jb08_openbook_g1
X_jb08 = X_true_aug(:,vec);
load starlink_550km_msis_openbook_g1
X_msis = X_true_aug(:,vec);
doy_g5 = doy;
err_msis = X_msis - X_hasdm;
err_jb08 = X_jb08 - X_hasdm;
for ii = 1:numel(vec)
    rhat = X_hasdm(1:3,ii)/norm(X_hasdm(1:3,ii));
vhat = X_hasdm(4:6,ii)/norm(X_hasdm(4:6,ii));
    Xhat = rhat;
Zhat = cross(rhat,vhat);
Yhat = cross(Zhat,Xhat);
Yhat = Yhat/norm(Yhat);
ECI2RIC = [Xhat';Yhat';Zhat'];
err_msis_ric_550_g1(:,ii) = ECI2RIC*err_msis(1:3,ii);
err_jb08_ric_550_g1(:,ii) = ECI2RIC*err_jb08(1:3,ii);

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

del_M = (M_hasdm(ii) - M_msis(ii));
del_s_msis_g1_550(ii) = del_M*v_hasdm/n_hasdm;

del_M = (M_hasdm(ii) - M_jb08(ii));
del_s_jb08_g1_550(ii) = del_M*v_hasdm/n_hasdm;
end

%%
fig = figure;
subplot(3,1,1)
plot(time_prop_utc/86400,del_s_msis_250/1e3,'-k','LineWidth',2)
hold on
plot(time_prop_utc/86400,del_s_jb08_250/1e3,'-r','LineWidth',2)
plot(time_prop_utc/86400,del_s_msis_g1_250/1e3,'--k','LineWidth',2)
plot(time_prop_utc/86400,del_s_jb08_g1_250/1e3,'--r','LineWidth',2)
grid on
legend('NRLMSISE-00','JB2008')
title('In-track error w.r.t HASDM (km) - 250 km');
set(gca,'FontSize',14)
yticks([-300 -200 -100 0 100 200 300])

subplot(3,1,2)
plot(time_prop_utc/86400,del_s_msis/1e3,'-k','LineWidth',2)
hold on
plot(time_prop_utc/86400,del_s_jb08/1e3,'-r','LineWidth',2)
plot(time_prop_utc/86400,del_s_msis_g1/1e3,'--k','LineWidth',2)
plot(time_prop_utc/86400,del_s_jb08_g1/1e3,'--r','LineWidth',2)
grid on
title('In-track error w.r.t HASDM (km) - 350 km');
set(gca,'FontSize',14)
yticks([-50 -25 0 25 50])
subplot(3,1,3)
plot(time_prop_utc/86400,del_s_msis_550/1e3,'-k','LineWidth',2)
hold on
plot(time_prop_utc/86400,del_s_jb08_550/1e3,'r','LineWidth',2)
plot(time_prop_utc/86400,del_s_msis_g1_550/1e3,'--k','LineWidth',2)
plot(time_prop_utc/86400,del_s_jb08_g1_550/1e3,'--r','LineWidth',2)
grid on
title('In-track error w.r.t HASDM (km) - 550 km');
set(gca,'FontSize',14)
yticks([-2 -1 0 1 2 3])
ylim([-2 3])
% han=axes(fig,'visible','off'); 
% han.Title.Visible='on';
% han.XLabel.Visible='on';
% han.YLabel.Visible='on';
% ylabel(han,'In-track error (km)');
xlabel('Days');
set(gca,'FontSize',14)

%% PLots
load DeltaV_dragCancellation_G1
delv1 = DelV_mat;
load DeltaV_dragCancellation_G5
delv2 = DelV_mat;
figure(1)
subplot(2,1,1)
plot(Hmat,delv1, 'k.', 'MarkerSize', 15)
hold on
load DeltaV_dragCancellation_G5
plot(Hmat,delv2, 'r.', 'MarkerSize', 15)
title('Delta v requirement over one day (m/s)')
legend('G1','G5')
grid on
set(gca,'FontSize',16)
subplot(2,1,2)
plot(Hmat,delv2./delv1, 'k.', 'MarkerSize', 15)
grid on
title('Ratio of required delta v (G5/G1)')
xlabel('Altitude (km)')
set(gca,'FontSize',16)
%%
load G5_dragseries
a_dragG1 = a_drag;
load G1_dragseries
a_dragG5 = a_drag;
figure(1)
plot(time_prop_utc_ekf/86400,a_dragG1*300*1000, 'k.')
hold on
load DeltaV_dragCancellation_G5
plot(time_prop_utc_ekf/86400,a_dragG5*300*1000, 'r.')
title('Drag force (mN)')
legend('G1','G5')
grid on
xlabel('Days')
set(gca,'FontSize',16)