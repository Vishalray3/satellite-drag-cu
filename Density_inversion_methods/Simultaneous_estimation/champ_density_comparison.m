%% Get density output from Champ 2006
clc
clearvars
Constants_ISS
load('champ_den_20061214')
rho_inputs = parameters;
epoch = -13;

r_circ = a_sma;
v_circ = sqrt(mu_e/r_circ);
n_mean = v_circ/r_circ;
T_orb = 2*pi/n_mean;

for ii = 1:numel(champ_mat(:,1))
    rho_inputs.altitude  = champ_mat(ii, 2)*1e3; % m
    rho_inputs.latitude  = champ_mat(ii, 3);
    rho_inputs.longitude = champ_mat(ii, 4);
    rho_inputs.time_jd   = jdutc(ii);
    rho_inputs.UTsec     = utc_time(ii);
    rho_inputs.doy       = doy_vec(ii);
    rho_inputs.sun_pos   = sun_pos(:,ii);
    rho_inputs.X_state   = x_eci(:,ii);
    rho(ii)                  = density_func(rho_inputs);
    
end

r_ind = round(T_orb/10);
r_alt = champ_mat(:, 2);
ii = 1;
time_prop_utc = time_prop_utc(1:end-1);
for kk=1:r_ind:numel(r_alt)-r_ind
    r_alt_avg(ii) = mean(r_alt(kk:kk+r_ind));
    r_alt_min(ii) = min(r_alt(kk:kk+r_ind));
    rho_avg(ii) = mean(rho(kk:kk+r_ind));
    rho_champ_avg(ii) = mean(champ_mat(kk:kk+r_ind,5));
    ii = ii+1;
end
time_avg = time_prop_utc(r_ind:r_ind:end);

ss = mod(time_prop_utc, 60);
min = mod(floor(time_prop_utc/60),60);
hr = mod(floor(time_prop_utc/3600),24);
dy = day_mon+floor(time_prop_utc/86400);
mn = mon*ones(1,numel(dy));
yr = yyyy*ones(1,numel(dy));

t_iss = datetime(yr,mn,dy,hr,min,ss);
t_issav = t_iss(r_ind:r_ind:end);
save('champ_hasdm_2006','rho_avg','t_iss','t_issav','rho_champ_avg', 'rho')
%%
time_vec = [0:10:259190];
figure(1)
load 'champ_msis_2006'
plot(t_issav,rho_avg,'LineWidth',2)
hold on
load 'champ_jb08_2006'
plot(t_issav,rho_avg,'LineWidth',2)
load 'champ_hasdm_2006'
plot(t_issav,rho_avg,'LineWidth',2)
plot(t_issav,rho_champ_avg,'k','LineWidth',4)
grid on
legend('NRLMSISE-00','JB2008','HASDM','CHAMP')
title('Orbit-averaged density (kg/m^3)');
set(gca,'FontSize',16)

figure(2)
load 'champ_hasdm_2006'
% plot(t_iss,champ_mat(:,5),'LineWidth',1)
hold on
plot(t_iss,rho,'LineWidth',1)
load 'champ_jb08_2006'
plot(t_iss,rho,'LineWidth',1)
% load 'champ_msis_2006'
% plot(t_iss,rho,'LineWidth',1)
grid on
legend('CHAMP','HASDM','JB2008','NRLMSISE-00')
title('Density (kg/m^3)');
set(gca,'FontSize',16)