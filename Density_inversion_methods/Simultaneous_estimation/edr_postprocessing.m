%% EDR Post processing
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Density inversion methods/Simultaneous estimation/error analysis/EDR')
clc
clearvars
alt_array = [3]*100;
t_array = [10, 60*[1:150]]; %[60, 60*[10:10:90], 60*[91:5:120]];
sample_time = 1;
for alt_index = 1:numel(alt_array)
    [rho_mean1,rho_rms1,rho_mean2,rho_rms2,T_orb, Ind] = edr_densities(alt_index,sample_time,t_array,alt_array);
    rho_est_mean_avg1(alt_index,:) = rho_mean1;
    rho_est_rms_avg1(alt_index,:) = rho_rms1;
    rho_est_mean_avg2(alt_index,:) = rho_mean2;
    rho_est_rms_avg2(alt_index,:) = rho_rms2;
    Torb_mat(alt_index) = T_orb;
    torb_ind(alt_index) = Ind;
     
end
save('EDR_10s_sampling_1cmNoise')
function [rho_mean1,rho_rms1,rho_mean2,rho_rms2,T_orb, Ind] = edr_densities(alt_index,sample_time,t_array,alt_array)
    load(strcat('EDR',num2str(alt_array(alt_index)), 'km_1s_sampling_10cmnoise'))
      
    for t_index = 1:numel(t_array)
        clear rho_true_eff  g_con_integ  g_srp_integ  g_erp_integ  g_noncon_integ  g_dragpar_integ  g_drag_integ  delta_v  delta_ener  rho_eff  time_rho rho_true_inst
    delta_t = t_array(t_index);
    kk = 1;
    for nn_ind = 1:sample_time:numel(time_interp)-delta_t
        x_state = X_interp(1:3,nn_ind:nn_ind+delta_t-1);
        v0 = norm(X_interp(4:6,nn_ind));
        vf = norm(X_interp(4:6,nn_ind+delta_t-1));
        g_con = a_grav(:,nn_ind:nn_ind+delta_t-1) + a_sun(:,nn_ind:nn_ind+delta_t-1) + a_moon(:,nn_ind:nn_ind+delta_t-1);
        g_srp = a_srp(:,nn_ind:nn_ind+delta_t-1);
        g_erp = a_earthrad(:,nn_ind:nn_ind+delta_t-1);
        g_dragpar = drag_par(:,nn_ind:nn_ind+delta_t-1);
        g_drag = a_drag(:,nn_ind:nn_ind+delta_t-1);
        
        rho_true_eff(kk) = mean(vecnorm(a_drag(:,nn_ind:nn_ind+delta_t-1),2,1))/mean(vecnorm(drag_par(:,nn_ind:nn_ind+delta_t-1),2,1));
        rho_true_inst(kk) = rho(nn_ind+round((delta_t-1)/2));
        
        g_con_integ(kk) = -(trapz(x_state(1,:),g_con(1,:))+ trapz(x_state(2,:),g_con(2,:))+ trapz(x_state(3,:),g_con(3,:)));
        g_srp_integ(kk) = trapz(x_state(1,:),g_srp(1,:))+ trapz(x_state(2,:),g_srp(2,:))+ trapz(x_state(3,:),g_srp(3,:));
        g_erp_integ(kk) = trapz(x_state(1,:),g_erp(1,:))+ trapz(x_state(2,:),g_erp(2,:))+ trapz(x_state(3,:),g_erp(3,:));
        g_noncon_integ(kk) = g_srp_integ(kk) + g_erp_integ(kk);
        
        g_dragpar_integ(kk) = trapz(x_state(1,:),g_dragpar(1,:))+ trapz(x_state(2,:),g_dragpar(2,:))+ trapz(x_state(3,:),g_dragpar(3,:));
        g_drag_integ(kk) = trapz(x_state(1,:),g_drag(1,:))+ trapz(x_state(2,:),g_drag(2,:))+ trapz(x_state(3,:),g_drag(3,:));
        
        delta_v(kk) = (vf^2-v0^2)/2;
        delta_ener(kk) = (vf^2-v0^2)/2 + g_con_integ(kk);
        
        rho_eff(kk) = (delta_ener(kk) - g_noncon_integ(kk))/g_dragpar_integ(kk);
        time_rho(kk) = time_interp(nn_ind);
        kk = kk+1;
    end
    
    rho_est_error1 = (rho_true_eff-rho_eff)./rho_true_eff*100;
    rho_est_error2 = (rho_true_inst-rho_eff)./rho_true_inst*100;
    rho_mean1(t_index) = mean(rho_est_error1);
    rho_rms1(t_index)  = rms(rho_est_error1);
    rho_mean2(t_index) = mean(rho_est_error2);
    rho_rms2(t_index)  = rms(rho_est_error2)
    [Mel, Ind] = min(abs(t_array-T_orb));
    end
end