function rho = density_jb08(epoch, doy,year, t, altitude, latitude,flattening, X_state, sun_pos, SOLdata, DTCdata, ind_sol_orig, ind_mag_orig)
rho = zeros(1,numel(t));
for kk = 1:numel(t)
UTsec = mod(epoch+t(kk),86400);
n_days = floor((epoch+t(kk))/86400);
dday = doy + n_days;
    
SUN = [0,0];
SAT = [0,0,0];


UThour  = UTsec/3600;
hour = floor(UThour);
minute = floor(UTsec/60) - hour*60;
sec = UTsec - hour*3600 - minute*60;
[month, day_mon, ~,~, ~] = days2mdh (year, dday);
MJD = Mjday(year,month,day_mon,hour,minute,sec);

ra_Sun  = atan2(sun_pos(2,kk), sun_pos(1,kk));
dec_Sun = atan2(sun_pos(3,kk), sqrt(sun_pos(1,kk)^2+sun_pos(2,kk)^2));
SUN(1)  = ra_Sun;
SUN(2)  = dec_Sun;


SAT(1) = atan2(X_state(2,kk), X_state(1,kk));                   % right-ascension
SAT(2) = geocentricLatitude(latitude(kk)*pi/180,flattening,'radians');         % F = flattening    latitude
SAT(3) = altitude(kk)/1000;                                    % altitude in km


ind_sol = ind_sol_orig+dday-2;                  % one day lag

SOL = SOLdata(:,ind_sol);
F10 = SOL(4);
F10B = SOL(5);
S10 = SOL(6);
S10B = SOL(7);

% USE 2 DAY LAG FOR M10 FOR JB2008
SOL = SOLdata(:,ind_sol-1);
XM10 = SOL(8);
XM10B = SOL(9);

% USE 5 DAY LAG FOR Y10 FOR JB2008
SOL = SOLdata(:,ind_sol-4);
Y10 = SOL(10);
Y10B = SOL(11);

ind_mag = ind_mag_orig+dday-1;
DTC = DTCdata(:,ind_mag);
ii = floor(hour)+3;
DSTDTC = DTC(ii);

[Temp,rho(kk),M,rho_oxa] = JB2008(MJD,SUN,SAT,F10,F10B,S10,S10B,XM10,XM10B,Y10,Y10B,DSTDTC);
T = Temp(2);
end