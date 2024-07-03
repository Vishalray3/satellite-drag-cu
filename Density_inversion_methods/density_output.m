%% Function that calculates continuous or discrete density
function [rho, rhodiff, M, rho_oxa, T] = density_output(time_jd, X_state, sun_pos, eps, parameters)
parameters.X_state = X_state;

if strcmp(parameters.flag_rho,'Exp')
    [rho,M,rho_oxa,T] = density_func(parameters);
    rhodiff = -1/parameters.H_scale*rho;
else
    TEI = parameters.TEI;
    
    X_ecef = TEI*X_state(1:3);
    X_lla = ecef2lla(X_ecef');
    latitude = X_lla(1);
    longitude = X_lla(2);
    altitude = X_lla(3);
    
    [Year,doy,h,m,UTsec] = JDtoGREGORIAN_vector(time_jd/86400);
    
    parameters.UTsec = UTsec;
    parameters.doy = doy;
    parameters.altitude = altitude;
    parameters.latitude = latitude;
    parameters.sun_pos = sun_pos;
    parameters.time_jd = time_jd;
    
    parameters.longitude = longitude;
    [rho, M, rho_oxa, T] = density_func(parameters);
    
    altp = altitude + eps;
    altn = altitude - eps;
    parameters.altitude = altp;
    rhop = density_func(parameters);
    parameters.altitude = altn;
    rhon = density_func(parameters);
    rhodiff = (rhop-rhon)/(2*eps);
end

end