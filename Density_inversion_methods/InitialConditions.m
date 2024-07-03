%% Initial conditions
%% Spire-like high altitude
% Hp = 510*1e3;           % perigee altitude
% Ha = 520*1e3;
% % a_sma= Re+H;          % only for circular
% r_p = Re+Hp;            % perigee radius
% r_a = Re+Ha;
% e = (r_a-r_p)/(r_a+r_p);   % eccentricity
% a_sma = r_p/(1-e);        % semi major axis from perigee altitude
% inc = 65;                  %%%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 97.8
% raan = 0;
% w_arg = 40;         % random
% true_ano = 0;      % random
% u_arg = 0;          % argument of latitude
% str_special = 'NO'; % CI: Circular inclined ; CE: circular equatorial; EE: elliptical equatorial; NO: None

%% CHAMP orbit on 29 Oct 2003
% Hp = 385.354*1e3;           % perigee altitude
% Ha = 397.860*1e3;
% % a_sma= Re+H;          % only for circular
% r_p = Re+Hp;            % perigee radius
% r_a = Re+Ha;
% e = (r_a-r_p)/(r_a+r_p);   % eccentricity
% a_sma = r_p/(1-e);        % semi major axis from perigee altitude
% inc = 87.266;                  %%%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 97.8
% raan = 51.6873;
% w_arg = 87.0484;         % random
% true_ano = 168.5641;      % random
% u_arg = -104.3875;          % argument of latitude
% str_special = 'NO'; % CI: Circular inclined ; CE: circular equatorial; EE: elliptical equatorial; NO: None

%% GRACE A orbit on 29 Oct 2003
% Hp = 460.823*1e3;           % perigee altitude
% Ha = 497.692*1e3;
% % a_sma= Re+H;          % only for circular
% r_p = Re+Hp;            % perigee radius
% r_a = Re+Ha;
% e = (r_a-r_p)/(r_a+r_p);   % eccentricity
% a_sma = r_p/(1-e);        % semi major axis from perigee altitude
% inc = 88.987;                  %%%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 97.8
% raan = -83.773;
% w_arg = 139.985;         % random
% true_ano = -178.318;      % random
% u_arg = -38.334;          % argument of latitude
% str_special = 'NO'; % CI: Circular inclined ; CE: circular equatorial; EE: elliptical equatorial; NO: None



%% Spire-like low altitude
Hp = Hp_ind*1e3;           % perigee altitude
Ha = Ha_ind*1e3;
% a_sma= Re+H;          % only for circular
r_p = Re+Hp;            % perigee radius
r_a = Re+Ha;
e = (r_a-r_p)/(r_a+r_p);   % eccentricity
a_sma = r_p/(1-e);        % semi major axis from perigee altitude
inc = 65;                  %%%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 97.8
raan = 0;
w_arg = 40;         % random
true_ano = 0;      % random
u_arg = 0;          % argument of latitude
str_special = 'NO'; % CI: Circular inclined ; CE: circular equatorial; EE: elliptical equatorial; NO: None



%% Starlink simulations
% Hp = 250*1e3;           % perigee altitude
% Ha = 250*1e3;
% % a_sma= Re+H;          % only for circular
% r_p = Re+Hp;            % perigee radius
% r_a = Re+Ha;
% e = (r_a-r_p)/(r_a+r_p);   % eccentricity
% a_sma = r_p/(1-e);        % semi major axis from perigee altitude
% inc = 0;                  %%%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 97.8
% raan = 0;
% w_arg = 40;         % random
% true_ano = 0;      % random
% u_arg = 0;          % argument of latitude
% str_special = 'CE'; % CI: Circular inclined ; CE: circular equatorial; EE: elliptical equatorial; NO: None

%% ISS
% inc = 51.6359;
% raan = 340.4161;
% e = 0.0023729;
% w_arg = 333.0463; 
% m_ano = 140.7139;
% mean_mot = 15.75851208461531*2*pi/86400;
% a_sma = (mu_e/mean_mot^2)^(1/3);
% str_special = 'NO';
% true_ano = m_ano;      % random
% u_arg = 0;          % argument of latitude