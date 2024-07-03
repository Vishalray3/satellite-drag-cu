%% Sample code to calculate drag-coefficient of a satellite with flat sufaces
clc
clearvars
%% User specified inputs
Alpha = 1;
Kl = 5e6;
frac_ox = 0.9;   
%% Satellite specific inputs
Tw = 300;                                           % Temperature of the satellite wall
% plate normal vectors
area_vec(:,1) = [1;0;0];
area_vec(:,2) = [-1;0;0];
area_vec(:,3) = [0;1;0];
area_vec(:,4) = [0;-1;0];
area_vec(:,5) = [0;0;1];
area_vec(:,6) = [0;0;-1];
area_vec(:,7) = sqrt(2)/2*[ 1;-1; 0];               % front solar panel half
area_vec(:,8) = sqrt(2)/2*[ 1;-1; 0];               % front solar panel half
area_vec(:,9) = sqrt(2)/2*[-1; 1; 0];               % Aft solar panel

% plate areas (m2)
Area_plates = [0.1053314*0.3388226, 0.1053314*0.3388226, 0.105186*0.3388226 + 3*0.005334*0.074549, 0.105186*0.3388226 + 3*0.005334*0.074549, ...
    0.105186*0.1053314, 0.105186*0.1053314, 0.150*0.3222, 0.150*0.3222, 2*0.150*0.3222];
Ar = 1;                                             % Reference area = 1

M_s = [26.981538, 60.0843, 60.0843, 26.981538, 26.981538, 26.981538, 26.981538, 60.0843, 60.0843];             % satellite surface material molar masses g/mol

% Quantities that'll utilize functions in Geodyn (attitude information, frame conversion, orbit propagation)
% Relative velocity vector in body frame  %% what coordinate frames are
% used in Geodyn
u = [1;0;0];
% Satellite relative speed
Vi = 7000;
%% Atmospheric inputs
atm_mass = [4.002602,15.9994,28.0134,31.9988,39.948,1.0079,14.0067, 15.9994]; % amu's or molar mass:g/mol; He, O, N2, O2, Ar, H, N, O

% Will need to call a density model in Geodyn to calculate the mean molar mass and
% ambient temperature
altitude = 300e3; latitude = 0; longitude = 0; year = 2018; day = 30; UTsec = 0; F10 = 75; F10d = 80; ap = [4,4,4,4,4,4,4]; 
[Temp,n_den] = atmosnrlmsise00(altitude, latitude, longitude, year, day, UTsec, F10, F10d, ap,'Oxygen');

mass_sum = n_den(:,1)*atm_mass(1) + n_den(:,2)*atm_mass(2) + n_den(:,3)*atm_mass(3) + n_den(:,4)*atm_mass(4) + n_den(:,5)*atm_mass(5) ...
          + n_den(:,7)*atm_mass(6) + n_den(:,8)*atm_mass(7) + n_den(:,9)*atm_mass(2);
num = sum(n_den,2);
M = mass_sum./num;           % Mean molar mass
n_o     = n_den(:,2) + n_den(:,9);    % density of atomic oxygen
T = Temp(2);

%% Calculate Cd

Cd_total = cd_dria_geodyn(u, area_vec, Tw, Vi, Area_plates, Ar, M_s, Alpha, M, T, n_o, Kl, frac_ox);