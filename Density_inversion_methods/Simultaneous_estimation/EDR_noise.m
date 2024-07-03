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
sigma_pos = sigma_pos_new;
sigma_vel = sigma_vel_new;
for ii = 1:numel(time_prop)
    Pos_noise(:,ii) =[sigma_pos*randn;sigma_pos*randn;sigma_pos*randn];
    Vel_noise(:,ii) = [sigma_vel*randn;sigma_vel*randn;sigma_vel*randn];
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
end
drag_par = a_drag./rho;

dragmod = mean(a_drag);
srpmod = mean(a_srp);
erpmod = mean(a_earthrad);
% save(truth_model)
% save(strcat('EDR', num2str(Halt_ind), 'km_1s_sampling_10cmnoise'))
%% Integrate the accelerations to calculate the energy decay
% clc
% clearvars
% load Truth800km_solarmin_1s
% time_interp = time_prop_utc;
% X_interp = X_true_aug(1:6,:);
% % clear rho_true_eff  g_con_integ  g_srp_integ  g_erp_integ  g_noncon_integ  g_dragpar_integ  g_drag_integ  delta_v  delta_ener  rho_eff  time_rho
arclen = [10, 60:60*60:120*60];
sample_time = 10;
for mm = 1:numel(arclen)
    mm
    delta_t = arclen(mm); %10; %round(T_orb);
    kk = 1;
    for ii = 1:sample_time:numel(time_interp)-delta_t
        x_state = X_interp(1:3,ii:ii+delta_t-1);
        v_state = X_interp(4:6,ii:ii+delta_t-1);
        v0 = norm(X_interp(4:6,ii));
        vf = norm(X_interp(4:6,ii+delta_t-1));
        g_con = a_grav(:,ii:ii+delta_t-1) + a_sun(:,ii:ii+delta_t-1) + a_moon(:,ii:ii+delta_t-1);
        g_srp = a_srp(:,ii:ii+delta_t-1);
        g_erp = a_earthrad(:,ii:ii+delta_t-1);
        g_dragpar = drag_par(:,ii:ii+delta_t-1);
        g_drag = a_drag(:,ii:ii+delta_t-1);
        
        rho_true_eff(kk) = mean(vecnorm(a_drag(:,ii:ii+delta_t-1),2,1))/mean(vecnorm(drag_par(:,ii:ii+delta_t-1),2,1));
        rho_true_inst(kk) = rho(ii+round((delta_t-1)/2));
        
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
        rho_true_eff_vec(ii:ii+delta_t-1) = rho_true_eff(kk);
        time_ind = ii + round(delta_t/2);
        time_rho(kk) = time_interp(time_ind);
        time_eff_vec(ii:ii+delta_t-1) = time_interp(ii:ii+delta_t-1);
        %% Calculate numerical integration error
        a_force = g_con + g_srp + g_erp;
        d2yd2x = diff(a_force.*v_state,2,2);
        max_d2yd2x = max(abs(d2yd2x), [], 2).*(delta_t).^3;
        num_points = delta_t;
        error_rhs_abs(kk) = abs(sum(abs(max_d2yd2x))/(12*num_points^2));
        error_rhs = abs(error_rhs_abs(kk) / (g_con_integ(kk) + g_noncon_integ(kk)));
        
        d2yd2x = diff(g_dragpar.*v_state,2,2);
        max_d2yd2x = max(abs(d2yd2x), [], 2).*(delta_t).^3;
        error_lhs = abs(sum(abs(max_d2yd2x))/(12*num_points^2) / g_dragpar_integ(kk));
        
%         error_rho_rel(kk) = (error_lhs + error_rhs)*100;
%         error_rho_abs(kk) = (error_lhs + error_rhs)*abs((delta_ener(kk) - g_noncon_integ(kk))/g_dragpar_integ(kk));
%         error_rho_rel(kk) = error_rho_abs(kk)/rho_true_inst(kk)*100;

        delta_ener_drag(kk) = g_dragpar_integ(kk)*rho_true_eff(kk);
        kk = kk+1;
        
    end
    
    rho_est_error = (rho_true_inst-rho_eff)./rho_true_inst*100;
    rho_est_error_avg = (rho_true_eff-rho_eff)./rho_true_eff*100;
    rho_est_mean_inst(mm) = mean(rho_est_error);
    rho_est_rms_inst(mm)  = rms(rho_est_error);
    rho_est_mean_avg(mm) = mean(rho_est_error_avg);
    rho_est_rms_avg(mm)  = rms(rho_est_error_avg);
    ener_drag_mean(mm) = mean(delta_ener_drag);
    rhs_abs_mean(mm) = mean(error_rhs_abs);
    
    clear rho_true_eff g_con_integ g_srp_integ g_erp_integ g_noncon_integ g_dragpar_integ g_drag_integ delta_ener rho_eff time_rho rho_true_inst error_rho_rel 
end

[rms_min, ind_min] = min(rho_est_rms_inst);
rho_est_rms = rms_min;
rho_est_mean = rho_est_mean_inst(ind_min);
P_signal = arclen(ind_min);
rms_min = min(rho_est_rms_avg);
ind_min = find(rho_est_rms_avg < 1.1*rms_min);
rho_rms_avg_min = rho_est_rms_avg(ind_min(1));
rho_mean_avg_min = rho_est_mean_avg(ind_min(1));
arc_min = arclen(ind_min(1));
ind_10 = find(rho_est_rms_avg < 10);
if ~isempty(ind_10)
    arc_10 = arclen(ind_10(1));
else
    arc_10 = NaN;
end
ind_5 = find(rho_est_rms_avg < 5);
if ~isempty(ind_5)
    arc_5 = arclen(ind_5(1));
else
    arc_5 = NaN;
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
