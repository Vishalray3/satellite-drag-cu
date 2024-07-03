%% Get Ap and F10.7 from Eric's tables

load formatted_geomag_1960_2019

year = 2006;

index_y = (geomag.yyyy == year);

Ap_daily = geomag.Ap_daily(index_y)';
F10_total = geomag.F107_daily(index_y)';
ap = geomag.ap(:,index_y);
Ap_total = ap(:);

save(strcat('NRLMSISE_', num2str(year)), 'Ap_daily','F10_total','Ap_total')