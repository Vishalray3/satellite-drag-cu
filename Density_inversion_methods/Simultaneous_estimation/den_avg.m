function [rho_est_mean_avg, rho_est_rms_avg] = den_avg(rho_true, rho_est, arclen) 

N = numel(rho_true);
for mm = 1:numel(arclen)
    delta_t = arclen(mm)/10; %10; %round(T_orb);
    kk = 1;
    for ii = 1:N-delta_t       
        rho_true_eff(kk) = mean(rho_true(ii:ii+delta_t-1));
        rho_eff(kk) = mean(rho_est(ii:ii+delta_t-1));
        kk = kk+1;
    end
    
    rho_est_error_avg = (rho_true_eff-rho_eff)./rho_true_eff*100;
    rho_est_mean_avg(mm) = mean(rho_est_error_avg);
    rho_est_rms_avg(mm)  = rms(rho_est_error_avg);
    clear rho_true_eff rho_eff
end
