%% Read CHAMP data
clc
clearvars
daysofyear = [348:350];
n_days = numel(daysofyear);
champ_mat = [];
utc_time = [];
jdutc = [];
yy = [2006,2006,2006];
month = [12,12,12];
dayofmon = [14, 15, 16];
doy_vec = [];
for jj = 1:numel(daysofyear)
den_mat = readmatrix(strcat('CHAMP_Density_06_', num2str(daysofyear(jj)),'_v1'));
% den_mat = readmatrix('graceA_Density_03_302_v1');
gps_time = den_mat(:,1);
tai_time = gps_time + 19;
utc_time1 = tai_time - 32;   %% for 2003
jdutc1 = 86400*GREGORIANtoJD_vector(yy(jj),month(jj),dayofmon(jj)) +utc_time1;
jd_ref = 86400*GREGORIANtoJD_vector(1950,1,0) + 3600*0 + 60*0 + 0;
njd_utc = jdutc1 - jd_ref;

champ_mat_temp(:,1) = njd_utc;
champ_mat_temp(:,2) = den_mat(:,2); % altitude in km
champ_mat_temp(:,3) = den_mat(:,3); % latitude
champ_mat_temp(:,4) = den_mat(:,4); % longitude
champ_mat_temp(:,5) = den_mat(:,16); % density
champ_mat = [champ_mat;champ_mat_temp];
utc_time = [utc_time; utc_time1];
jdutc  = [jdutc; jdutc1];
doy_vec = [doy_vec; daysofyear(jj)*ones(8640,1)];
end
%% Find position vector
load nut80.dat
eop = read_finals_data;
eop.dat = read_tai_utc_dat; % reads only the most current value in the local data file
eqeterms = 1;
% ecef position
ecef = lla2ecef([champ_mat(:,3), champ_mat(:,4), champ_mat(:,2)*1e3]);

% eci
for ii = 1:numel(jdutc)
[rot_ECI2ECEF,dut1,xp,yp] = time2rotmat(eop, jdutc(ii),[0;0;0;0;0;0], eqeterms, nut80);
x_eci(:,ii) = rot_ECI2ECEF'*ecef(ii,:)';
end

%% Gibbs method to find velocity
mu_e = 0.3986004415E+15;
Re = 6378136.3;
r1 = x_eci(:,400);
r2 = x_eci(:,555);
r3 = x_eci(:,700);

Z12 = cross(r1,r2);
Z23 = cross(r2,r3);
Z31 = cross(r3,r1);
alpha_12 = acosd(dot(r1/norm(r1), r2/norm(r2)));
alpha_23 = acosd(dot(r2/norm(r2), r3/norm(r3)));

N_vec = norm(r1)*Z23 + norm(r2)*Z31 + norm(r3)*Z12;
D_vec = Z12+Z23+Z31;
S_vec = (norm(r2) - norm(r3))*r1 + (norm(r3) - norm(r1))*r2 + (norm(r1) - norm(r2))*r3;
B_vec = cross(D_vec,r2);
Lg = sqrt(mu_e/(norm(N_vec)*norm(D_vec)));
v2 = Lg/norm(r2)*B_vec + Lg*S_vec;

coe = rv2coe_E(r2,v2,mu_e);
a_sma = coe(1); e = coe(2); inc = coe(3); raan = coe(4); w_arg = coe(5); true_ano = coe(6); u_arg = coe(7);
H_p = a_sma*(1-e) - Re;
H_a = a_sma*(1+e) - Re;

% %% Run constants with HASDM after generating the truth and then run this section
% for ii = 1:numel(champ_mat(:,1))
%     if champ_mat(ii,4) < 0 
%             champ_mat(ii,4) = champ_mat(ii,4) + 360;
%     end
%     rho_hasdm(ii) = exp(parameters.F_hasdm(champ_mat(ii,2), champ_mat(ii,1)/86400, champ_mat(ii,4), champ_mat(ii,3)));
% end
