clc
clearvars
load('spire_sat_sp383_2018_11_7','Xsp3_eci','sod_sec_sp3')
X_eci_true = Xsp3_eci;
sod_true = sod_sec_sp3;

load('spire_sat83_2018_11_7_corrected','Xsp3_eci','sod_sec_sp3')

X_interp = interp1(sod_true, X_eci_true',sod_sec_sp3,'spline');
X_interp = X_interp';

err = X_interp - Xsp3_eci;