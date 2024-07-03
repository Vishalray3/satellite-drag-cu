%% EDR calculations in ECEF frame 
%% EDR method
% clc
% clearvars
% addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Density inversion methods/Simultaneous estimation/error analysis')
% addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Density inversion methods')
% run Main_simulStudy

% clear estimated_coeff
% parameters.roll = 0;
% parameters.pitch = 0;
% parameters.yaw = 0;
% estimated_coeff.Cr = 0;
% estimated_coeff.Cd = 0;
% estimated_coeff.rho_DMC = 0;
% estimated_coeff.CdTrue = 0;
% estimated_coeff.CdDiscrete = 0;
% estimated_coeff.CdDiscrete_est = 0;
% estimated_coeff.Cr_erp = 0;
% parameters.estimated_coeff = estimated_coeff;
%
% parameters.sun       = 0;
% parameters.moon      = 0;
% parameters.srp       = 0;
% parameters.drag      = 1;
% parameters.tides     = 0;
% parameters.relativity= 0;
% flag_tides.SolidEarthTides = 0;
% flag_tides.OceanTides = 0;
% parameters.earth_rad = 0;
% parameters.empirical = 0;
%
% parameters.flag_earthrad = 'Knocke';
% parameters.Cr_erp = 1;


%% Interpolate the ephemeris to one second interval

% run Main_simulStudy
clc
clearvars
load Truth800km_edr
sigma_pos = 0.01; %sigma_pos_new;
sigma_vel = 1e-5; %sigma_vel_new;
for ii = 1:numel(time_prop)
    Pos_noise(:,ii) = [normrnd(0,sigma_pos);normrnd(0,sigma_pos);normrnd(0,sigma_pos)];
    Vel_noise(:,ii) = [normrnd(0,sigma_vel);normrnd(0,sigma_vel);normrnd(0,sigma_vel)];
    GPS_state = X_true_aug(1:6,ii);
    reci(:,ii) = GPS_state(1:3) + Pos_noise(:,ii);
    veci(:,ii) = GPS_state(4:6) + Vel_noise(:,ii);
end
X_true_noise(1:6,:) = [reci;veci];
time_interp   = [time_prop_utc(1):time_prop_utc(end)];
X_interp(1,:) = interp1(time_prop_utc, X_true_noise(1,:), time_interp, 'spline');
X_interp(2,:) = interp1(time_prop_utc, X_true_noise(2,:), time_interp, 'spline');
X_interp(3,:) = interp1(time_prop_utc, X_true_noise(3,:), time_interp, 'spline');
X_interp(4,:) = interp1(time_prop_utc, X_true_noise(4,:), time_interp, 'spline');
X_interp(5,:) = interp1(time_prop_utc, X_true_noise(5,:), time_interp, 'spline');
X_interp(6,:) = interp1(time_prop_utc, X_true_noise(6,:), time_interp, 'spline');
time_interp_jdutc = [time_prop(1):time_prop(end)];

X_interp_aug = zeros(numel(X_true_noise(:,1)),numel(time_interp));
X_interp_aug(1:6,:) = X_interp;

% spec_ref = [0.6, 0.5, 0.5, 0.5, 0.8, 0.1, 0.7, 0.6, 0.9];
% diff_ref = [0.8, 0.7, 0.1, 0.2, 0.6, 0.9, 0.2, 0.1, 0.6];
% spec_ref_ir = [0.2, 0.05, 0.05, 0.2, 0.2, 0.2, 0.2, 0.05, 0.05];
% diff_ref_ir = [0.4, 0.3, 0.3, 0.4, 0.4, 0.4, 0.4, 0.3, 0.3];
% parameters.rho_diff = diff_ref;
% parameters.rho_spec = spec_ref;
% parameters.rho_diff_ir = diff_ref_ir;
% parameters.rho_spec_ir = spec_ref_ir;

% flag_earthrad = 'Knocke';
% parameters.flag_earthrad = flag_earthrad;
% parameters.roll = roll_ind;
% parameters.pitch = pitch_ind;
% parameters.yaw = pitch_ind;
% Calculate the different accelerations
for ii = 1:numel(time_interp)
    [~,Cd_est(ii),rho(ii), theta(ii), phi(ii),a_grav(:,ii),a_sun(:,ii),a_moon(:,ii),a_drag(:,ii),a_srp(:,ii),a_earthrad(:,ii)] = ...
        propagator_truth(time_interp(ii),X_interp_aug(:,ii),parameters);
    aeci(:,1) = a_sun(:,ii) + a_moon(:,ii) + a_srp(:,ii) + a_earthrad(:,ii);
    aeci(:,2) = a_drag(:,ii); 
    [rot_ECI2ECEF,dut1,xp,yp, recef(:,ii), vecef(:,ii), ~, omegae(ii)] = time2rotmat(eop, time_interp_jdutc(ii), X_interp_aug(:,ii),...
        parameters.eqeterms, parameters.nut80, aeci);
    a_ecef_int(:,ii) = rot_ECI2ECEF*aeci(:,1);
    drag_par(:,ii) = rot_ECI2ECEF*aeci(:,2)/rho(ii);   
end
Vecef = calculate_geopotential(parameters.deg_grav, parameters.ord_grav, Re, mu_e, recef, parameters.Cnm, parameters.Snm);
dragmod = mean(a_drag);
srpmod = mean(a_srp);
erpmod = mean(a_earthrad);
save(truth_model)
%% Integrate the accelerations to calculate the energy decay
clc
clearvars
load Truth800km_edr
% clear rho_true_eff  g_con_integ  g_srp_integ  g_erp_integ  g_noncon_integ  g_dragpar_integ  g_drag_integ  delta_v  delta_ener  rho_eff  time_rho
arclen = [10, 60, 10*60, [30:30:1200]*60];
sample_time = 10;
for mm = 1:numel(arclen)
delta_t = arclen(mm); %10; %round(T_orb);
kk = 1;
for ii = 1:sample_time:numel(time_interp)-delta_t
    x_state = recef(:,ii:ii+delta_t-1);
    v0 = norm(vecef(:,ii));
    vf = norm(vecef(:,ii+delta_t-1));
    g_con = a_grav(:,ii:ii+delta_t-1) + a_sun(:,ii:ii+delta_t-1) + a_moon(:,ii:ii+delta_t-1);
    g_noncon = a_ecef_int(:,ii:ii+delta_t-1);
    g_dragpar = drag_par(:,ii:ii+delta_t-1);
    
    rho_true_eff(kk) = mean(vecnorm(a_drag(:,ii:ii+delta_t-1),2,1))/mean(vecnorm(drag_par(:,ii:ii+delta_t-1),2,1));
    rho_true_inst(kk) = rho(ii+round((delta_t-1)/2));
    
    g_con_integ(kk) =  -(Vecef(ii+delta_t-1) - Vecef(ii)) - ...
        1/2*( (omegae(ii+delta_t-1)^2*(recef(1,ii+delta_t-1)^2 + recef(2,ii+delta_t-1)^2)) -  (omegae(ii)^2*(recef(1,ii)^2 + recef(2,ii)^2)) );
    g_noncon_integ(kk) = trapz(x_state(1,:),g_noncon(1,:))+ trapz(x_state(2,:),g_noncon(2,:))+ trapz(x_state(3,:),g_noncon(3,:));
    
    g_dragpar_integ(kk) = trapz(x_state(1,:),g_dragpar(1,:))+ trapz(x_state(2,:),g_dragpar(2,:))+ trapz(x_state(3,:),g_dragpar(3,:));
    %g_drag_integ(kk) = trapz(x_state(1,:),g_drag(1,:))+ trapz(x_state(2,:),g_drag(2,:))+ trapz(x_state(3,:),g_drag(3,:));
    
    delta_v(kk) = (vf^2-v0^2)/2;
    delta_ener(kk) = (vf^2-v0^2)/2 + g_con_integ(kk);
    
    rho_eff(kk) = (delta_ener(kk) - g_noncon_integ(kk))/g_dragpar_integ(kk);
    rho_eff_vec(ii:ii+delta_t-1) = rho_eff(kk);
    rho_true_eff_vec(ii:ii+delta_t-1) = rho_true_eff(kk);
    time_ind = ii + round(delta_t/2);
    time_rho(kk) = time_interp(time_ind);
    time_eff_vec(ii:ii+delta_t-1) = time_interp(ii:ii+delta_t-1);
    kk = kk+1;
end

rho_est_error = (rho_true_inst-rho_eff)./rho_true_inst*100;
rho_est_error_avg = (rho_true_eff-rho_eff)./rho_true_eff*100;
rho_est_mean(mm) = mean(rho_est_error);
rho_est_rms(mm)  = rms(rho_est_error);
rho_est_mean_avg(mm) = mean(rho_est_error_avg);
rho_est_rms_avg(mm)  = rms(rho_est_error_avg);
clear rho_true_eff g_con_integ g_srp_integ g_erp_integ g_noncon_integ g_dragpar_integ g_drag_integ delta_ener rho_eff time_rho rho_true_inst
end
% save('EDR_1s_sampling_10cmnoise')
%% PLots
% figure(1)
% subplot(2,1,1)
% % load('EDR_300km_10s')
% % plot(time_interp/3600,rho,'k','LineWidth',2)
% plot(time_rho/3600,rho_true_inst,'k','LineWidth',2)
% hold on
% plot(time_rho/3600,rho_eff,'r','LineWidth',2)
% grid on
% ylabel('Density ($kg/m^3$)','Interpreter','latex')
% title('Estimated density (300 km)')
% legend('Truth (HASDM)','Estimate')
% set(gca,'FontSize',16)

% subplot(2,1,2)
% % load('EDR_800km_Torb')
% plot(time_interp/3600,rho,'k','LineWidth',2)
% hold on
% plot(time_eff_vec/3600,rho_true_eff_vec,'b','LineWidth',2)
% plot(time_eff_vec/3600,rho_eff_vec,'--r','LineWidth',2)
% grid on
% ylabel('Density ($kg/m^3$)','Interpreter','latex')
% legend('Truth (HASDM)','Averaged truth','Estimate')
% title('Estimated density (800 km)')
% set(gca,'FontSize',16)
% xlabel('Time (hours)')


