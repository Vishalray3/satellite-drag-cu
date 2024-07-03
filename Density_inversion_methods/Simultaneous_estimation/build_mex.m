%% Build mex functions 
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods')
clc
clearvars
num = 1;
time_var = coder.typeof(num, 1);
X_state = coder.typeof(num, [6 1]);
eop.yy = coder.typeof(num, [1 Inf]);
eop.mm = coder.typeof(num, [1 Inf]);
eop.dd = coder.typeof(num, [1 Inf]);
eop.fmjd = coder.typeof(num, [1 Inf]);
eop.xp = coder.typeof(num, [1 Inf]);
eop.yp = coder.typeof(num, [1 Inf]);
eop.dut1 = coder.typeof(num, [1 Inf]);
eop.rlod = coder.typeof(num, [1 Inf]);
eop.dpsi = coder.typeof(num, [1 Inf]);
eop.deps = coder.typeof(num, [1 Inf]);
eop.dat = coder.typeof(num, 1);
eqeterms = coder.typeof(num, 1);
nut80 = coder.typeof(num, [106 10]);

codegen -report time2rotmat.m -args {eop, time_var,X_state, eqeterms, nut80 } %-test test_mex

C = coder.typeof(num, [Inf Inf]);
S = coder.typeof(num, [Inf Inf]);
time_et = coder.typeof(num, 1);
Mjd_UTC = coder.typeof(num, 1);
UT1_UTC = coder.typeof(num, 1);
TEI = coder.typeof(num, [3 3]);
pos_Moon = coder.typeof(num, [3 1]);
pos_Sun = coder.typeof(num, [3 1]);
r_ref = coder.typeof(num, 1);
gm = coder.typeof(num, 1);
mu_m = coder.typeof(num, 1);
mu_s = coder.typeof(num, 1);
x_pole = coder.typeof(num, 1);
y_pole = coder.typeof(num, 1);
flag_tides.SolidEarthTides = coder.typeof(num, 1);
flag_tides.OceanTides  = coder.typeof(num, 1);
coeff0  = coder.typeof(num, [48 7]);
coeff1  = coder.typeof(num, [21 7]);
coeff2  = coder.typeof(num, [2 6]);

codegen -report tides.m -args {C, S, time_et, Mjd_UTC, UT1_UTC, TEI, pos_Moon, pos_Sun, r_ref, gm, mu_m, mu_s, x_pole, y_pole, flag_tides, coeff0, coeff1, coeff2} -test test_mex
