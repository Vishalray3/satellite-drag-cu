function [rho_true_results, rho_nom_results, rho_est_results, alt_results] = den_avg_spire(rho_true, rho_est, rho_nom, arclen, time_step, time_ind, alt)
%% Inputs
% arclen = arc-length in sec
% time_step = time steps at which densities are estimated


delta_ind = round(arclen/time_step);

rho_true_results = [];
rho_nom_results = [];
rho_est_results = [];
alt_results = [];

for mm = 1:numel(time_ind)-1
    rho_true_arc = rho_true(time_ind(mm):time_ind(mm+1));
    rho_est_arc = rho_est(time_ind(mm):time_ind(mm+1));
    rho_nom_arc = rho_nom(time_ind(mm):time_ind(mm+1));
    alt_arc = alt(time_ind(mm):time_ind(mm+1));
    
%     if delta_ind < numel(rho_true_arc)
        rho_true_avg = movmean(rho_true_arc, delta_ind,'Endpoints','discard');
        rho_nom_avg = movmean(rho_nom_arc, delta_ind,'Endpoints','discard');
        rho_est_avg = movmean(rho_est_arc, delta_ind,'Endpoints','discard');
        alt_avg = movmean(alt_arc, delta_ind,'Endpoints','discard');
%     else
%         rho_true_avg = mean(rho_true_arc);
%         rho_nom_avg = mean(rho_nom_arc);
%         rho_est_avg = mean(rho_est_arc);
%         alt_avg = mean(alt_arc);
%     end
     
    rho_true_results = [rho_true_results, rho_true_avg];
    rho_nom_results  = [rho_nom_results, rho_nom_avg];
    rho_est_results  = [rho_est_results, rho_est_avg];
    alt_results = [alt_results, alt_avg];
    
end


