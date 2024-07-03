%% test mex
Constants_simSpire
[rot_ECI2ECEF,dut1,xp,yp] = time2rotmat(eop, time_prop(1),X_init, eqeterms, nut80);
deg_grav = 80;
ord_grav = 80;
parameters.deg_grav = deg_grav;
parameters.ord_grav = ord_grav;
parameters.Cnm = Cbar(1:deg_grav+1, 1:deg_grav+1);
parameters.Snm = Sbar(1:deg_grav+1, 1:deg_grav+1);

parameters.sun_pos = sun_pos;
parameters.moon_pos = moon_pos;
[Cnm, Snm] = tides(parameters.Cnm, parameters.Snm, time_prop(1), time_prop(1), dut1, rot_ECI2ECEF, parameters.moon_pos(:,1), parameters.sun_pos(:,1), ...
    parameters.Re, mu_e,parameters.mu_moon, parameters.mu_sun, xp, yp, parameters.flag_tides, parameters.coeff0, parameters.coeff1, parameters.coeff2);