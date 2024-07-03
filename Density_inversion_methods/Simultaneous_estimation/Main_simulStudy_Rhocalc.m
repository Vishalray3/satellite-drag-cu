
%% Flags
if strcmp(dataset_sat, 'spire')
    Constants_spire
    propagator_dmc_real = @propagator_dmc_real_spire;
elseif strcmp(dataset_sat, 'starlink')
    Constants_starlink
    propagator_dmc_real = @propagator_dmc_real_starlink;
end


%             run ekf_dmc_smoothing_real
parameters.sun       = 0;
parameters.moon      = 0;
parameters.tides     = 0;
parameters.relativity= 0;
flag_tides.SolidEarthTides = 0;
flag_tides.OceanTides = 0;
deg_grav = 0;
parameters.deg_grav = deg_grav;
parameters.ord_grav = deg_grav;
parameters.Cnm = Cbar(1:deg_grav+1, 1:deg_grav+1);
parameters.Snm = Sbar(1:deg_grav+1, 1:deg_grav+1);
%             parameters.vec_est = [1:9];

if strcmp(dataset_sat, 'spire')
    load(strcat(dir_data,'/spireResults_', flag_rho, '_',sat_ID,'_', date),'time_new','X_full_uni','parameters','time_prop_utc_ekf','jdutc_sec_data','rho_est', 'rho_nom')
    for ii = 1:numel(time_new)
        [~,~,~,rho_nom_new(ii),~,~, ~,rot_ECI2ECEF, ~, ~] = propagator_dmc_real(time_new(ii),X_full_uni(:,ii),parameters);
        X_ecef = rot_ECI2ECEF*X_full_uni(1:3,ii);
        X_lla = ecef2lla(X_ecef');
        latitude(ii) = X_lla(1);
        longitude(ii) = X_lla(2);
        alt_calc(ii) = X_lla(3);
    end
%     F_den_est = griddedInterpolant(time_prop_utc_ekf, log(rho_est));
%     F_den_nom = griddedInterpolant(time_prop_utc_ekf, log(rho_nom));
    jdutc_new = interp1(time_prop_utc_ekf, jdutc_sec_data, time_new, 'linear', 'extrap');
    save(strcat(dir_data,'/spireResultsRho_', flag_rho, '_',sat_ID,'_', date),'latitude','longitude','alt_calc','jdutc_new')
elseif strcmp(dataset_sat, 'starlink')
    load(strcat(dir_data,'/starlinkResults_', flag_rho, '_',sat_ID,'_', date),'time_new','X_full_uni','parameters','time_prop_utc','jdutc_sec')
    for ii = 1:numel(time_new)
        [~,~,rho_est(ii),rho_nom(ii),~,~, ~,rot_ECI2ECEF, ~, ~] = propagator_dmc_real(time_new(ii),X_full_uni(:,ii),parameters);
        X_ecef = rot_ECI2ECEF*X_full_uni(1:3,ii);
        X_lla = ecef2lla(X_ecef');
        latitude(ii) = X_lla(1);
        longitude(ii) = X_lla(2);
        alt_calc(ii) = X_lla(3);        
    end
%     F_den_est = griddedInterpolant(time_new, log(rho_est));
%     F_den_nom = griddedInterpolant(time_new, log(rho_nom));
    jdutc_new = interp1(time_prop_utc, jdutc_sec, time_new, 'linear', 'extrap');
    save(strcat(dir_data,'/starlinkResultsRho_', flag_rho, '_',sat_ID,'_', date),'latitude','longitude','alt_calc','jdutc_new')
end
