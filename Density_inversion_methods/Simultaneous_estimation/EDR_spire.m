%% EDR for Spire data


%% Interpolate the ephemeris to one second interval and fill data-gaps with Nan
X_true_noise(1:6,:) = [reci;veci];
time_step = median(diff(time_prop_utc_ekf));
t1 = time_prop_utc_ekf(1);     % starting time
t_count = round((time_prop_utc_ekf-t1)/time_step);
t_index = t_count + 1;            % index of these times in array of regular times
ti = t1 + time_step * t_count;     % time to nearest timeStep


xi = interp1(time_prop_utc_ekf, X_true_noise', ti, 'spline');
xi = xi';

% evenly spaced times for whole sequence
time_interp = linspace(ti(1), ti(end), t_index(end));
% output array of x values
X_interp = NaN(6, t_index(end));
X_interp(:,t_index) = xi;

X_interp_aug = zeros(numel(X_true_noise(:,1)),numel(time_interp));
X_interp_aug(1:6,:) = X_interp;

%% Calculate the different accelerations
for ii = 1:numel(time_interp)
    [~,Cd_est(ii),rho(ii), rho_nom(ii), ~, a_drag(:,ii), ~, ~, ~, ~, rot_ECI2ECEF, a_grav(:,ii), a_sun(:,ii), a_moon(:,ii), a_srp(:,ii), a_earthrad(:,ii)]...
        = propagator_dmc_real_spire(time_interp(ii),X_interp_aug(:,ii),parameters);
    X_ecef = rot_ECI2ECEF*X_interp_aug(1:3,ii);
    X_lla = ecef2lla(X_ecef');
    latitude(ii) = X_lla(1);
    longitude(ii) = X_lla(2);
    alt_calc(ii) = X_lla(3);
end
drag_par = a_drag./rho;

dragmod = mean(a_drag);
srpmod = mean(a_srp);
erpmod = mean(a_earthrad);

%% Hasdm interpolation
F_hasdm = hasdm_interpolant(n_days, doy, hasdm_mat);
jd_ref = 86400*GREGORIANtoJD_vector(1950,1,0) + 3600*0 + 60*0 + 0;
jdutc_interp = jdutc_sec_data(1) + time_interp;

njd = (jdutc_interp - jd_ref)/86400;
long_pos = longitude;
long_pos(long_pos<0) = long_pos(long_pos<0) + 360;
rho_hasdm_interp = exp(F_hasdm(alt_calc/1e3, njd,long_pos, latitude));

a_drag_hasdm = rho_hasdm_interp.*drag_par;
%% Integrate the accelerations to calculate the energy decay
% clc
% clearvars
% load Truth800km_edr
% clear rho_true_eff  g_con_integ  g_srp_integ  g_erp_integ  g_noncon_integ  g_dragpar_integ  g_drag_integ  delta_v  delta_ener  rho_eff  time_rho
arclen = [60*60];
sample_time = 10;
for mm = 1:numel(arclen)
    delta_t = arclen(mm); %10; %round(T_orb);
    kk = 1;
    for ii = 1:sample_time:numel(time_interp)-delta_t
        x_state = X_interp(1:3,ii:ii+delta_t-1);
        v0 = norm(X_interp(4:6,ii));
        vf = norm(X_interp(4:6,ii+delta_t-1));
        g_con = a_grav(:,ii:ii+delta_t-1) + a_sun(:,ii:ii+delta_t-1) + a_moon(:,ii:ii+delta_t-1);
        g_srp = a_srp(:,ii:ii+delta_t-1);
        g_erp = a_earthrad(:,ii:ii+delta_t-1);
        g_dragpar = drag_par(:,ii:ii+delta_t-1);
        g_drag = a_drag(:,ii:ii+delta_t-1);
        
        rho_nom_eff(kk) = mean(vecnorm(a_drag(:,ii:ii+delta_t-1),2,1))/mean(vecnorm(drag_par(:,ii:ii+delta_t-1),2,1));
        rho_nom_inst(kk) = rho(ii+round((delta_t-1)/2));
        rho_hasdm_inst(kk) = rho_hasdm_interp(ii+round((delta_t-1)/2));
        
        rho_hasdm_eff(kk) = mean(vecnorm(a_drag_hasdm(:,ii:ii+delta_t-1),2,1))/mean(vecnorm(drag_par(:,ii:ii+delta_t-1),2,1));
        
        g_con_integ(kk) = -(trapz(x_state(1,:),g_con(1,:))+ trapz(x_state(2,:),g_con(2,:))+ trapz(x_state(3,:),g_con(3,:)));
        g_srp_integ(kk) = trapz(x_state(1,:),g_srp(1,:))+ trapz(x_state(2,:),g_srp(2,:))+ trapz(x_state(3,:),g_srp(3,:));
        g_erp_integ(kk) = trapz(x_state(1,:),g_erp(1,:))+ trapz(x_state(2,:),g_erp(2,:))+ trapz(x_state(3,:),g_erp(3,:));
        g_noncon_integ(kk) = g_srp_integ(kk) + g_erp_integ(kk);
        
        g_dragpar_integ(kk) = trapz(x_state(1,:),g_dragpar(1,:))+ trapz(x_state(2,:),g_dragpar(2,:))+ trapz(x_state(3,:),g_dragpar(3,:));
        g_drag_integ(kk) = trapz(x_state(1,:),g_drag(1,:))+ trapz(x_state(2,:),g_drag(2,:))+ trapz(x_state(3,:),g_drag(3,:));
        
        delta_v(kk) = (vf^2-v0^2)/2;
        delta_ener(kk) = (vf^2-v0^2)/2 + g_con_integ(kk);
        
        rho_eff(kk) = (delta_ener(kk) - g_noncon_integ(kk))/g_dragpar_integ(kk);
        rho_eff_vec(ii:ii+delta_t-1) = rho_eff(kk);
        rho_nom_eff_vec(ii:ii+delta_t-1) = rho_nom_eff(kk);
        rho_hasdm_eff_vec(ii:ii+delta_t-1) = rho_hasdm_eff(kk);
        time_ind = ii + round(delta_t/2);
        time_rho(kk) = time_interp(time_ind);
        time_eff_vec(ii:ii+delta_t-1) = time_interp(ii:ii+delta_t-1);
        kk = kk+1;
    end
    
    rho_est_error = (rho_hasdm_inst-rho_eff)./rho_hasdm_inst*100;
    rho_est_error_avg = (rho_hasdm_eff-rho_eff)./rho_hasdm_eff*100;
    rho_nom_error = (rho_hasdm_inst-rho_nom_eff)./rho_hasdm_inst*100;
    rho_nom_error_avg = (rho_hasdm_eff-rho_nom_eff)./rho_hasdm_eff*100;
    rho_est_mean(mm) = nanmean(rho_est_error);
    rho_est_rms(mm)  = sqrt(nanmean(rho_est_error.*rho_est_error));
    rho_est_mean_avg(mm) = nanmean(rho_est_error_avg);
    rho_est_rms_avg(mm)  = sqrt(nanmean(rho_est_error_avg.*rho_est_error_avg));
    rho_nom_mean(mm) = nanmean(rho_nom_error);
    rho_nom_rms(mm)  = sqrt(nanmean(rho_nom_error.*rho_nom_error));
    rho_nom_mean_avg(mm) = nanmean(rho_nom_error_avg);
    rho_nom_rms_avg(mm)  = sqrt(nanmean(rho_nom_error_avg.*rho_nom_error_avg));
%     clear rho_true_eff g_con_integ g_srp_integ g_erp_integ g_noncon_integ g_dragpar_integ g_drag_integ delta_ener rho_eff time_rho rho_nom_eff rho_hasdm_eff
end

% [rms_min, ind_min] = min(rho_est_rms);
% rho_est_rms = rms_min;
% rho_est_mean = rho_est_mean(ind_min);
% P_signal = arclen(ind_min);
% rms_min = min(rho_est_rms_avg);
% ind_min = find(rho_est_rms_avg < 1.1*rms_min);
% rho_rms_avg_min = rho_est_rms_avg(ind_min(1));
% rho_mean_avg_min = rho_est_mean_avg(ind_min(1));
% arc_min = arclen(ind_min(1));
% ind_10 = find(rho_est_rms_avg < 10);
% if ~isempty(ind_10)
%     arc_10 = arclen(ind_10(1));
% else
%     arc_10 = NaN;
% end
% ind_20 = find(rho_est_rms_avg < 20);
% if ~isempty(ind_20)
%     arc_20 = arclen(ind_20(1));
% else
%     arc_20 = NaN;
% end
save(strcat('results_cdpanel_spire_satellite_data',data_pod, '_', sat_ID,'_',date_curr))
%% PLots
figure(2)
plot(time_rho/3600,rho_hasdm_eff,'k','LineWidth',2)
hold on
plot(time_rho/3600,rho_nom_eff,'r','LineWidth',2)
plot(time_rho/3600,rho_eff,'g','LineWidth',2)
grid on
ylabel('Density ($kg/m^3$)','Interpreter','latex')
xlabel('Time (hours)')
title('Spire-derived densities (60-minute arc-length)')
legend('HASDM','JB2008','EDR')
set(gca,'FontSize',14)

