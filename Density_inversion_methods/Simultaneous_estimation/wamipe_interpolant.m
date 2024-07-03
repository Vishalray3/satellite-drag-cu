%% Read WAM-IPE data and create density interpolant
function F_wamipe = wamipe_interpolant(filename, year_data, month_data, day_data)
load(filename)
[hh, mm, ss] = hms([0:minutes(10):minutes(1430)]);
yyyy = repmat(year_data, 1, numel(hh));
mon = repmat(month_data, 1, numel(hh));
dy = repmat(day_data, 1, numel(hh));
time_jd = 86400*GREGORIANtoJD_vector(yyyy,mon,dy) + 3600*hh + 60*mm + ss;
jd_ref = 86400*GREGORIANtoJD_vector(1950,1,0) + 3600*0 + 60*0 + 0;
vec_njd = (time_jd - jd_ref)/86400;

gridVecs = {double(vec_long),double(vec_lat),double(vec_alt),vec_njd'};
F_wamipe = griddedInterpolant(gridVecs,double(rho_wam));
end
