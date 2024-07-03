function Cd_total = cd_dria_geodyn(u, Ai, Tw, Vi, A, Ar, M_s, Alpha, M, T, n_oxa, Kl, frac_ox)
%% Function to calculate physics-based value of drag-coefficients using a panel based method
% Assumptions -
% Satellite composed of flat panels. Sums up the drag-coefficient for individual panels to calculate the total drag-coefficient

% Model used -  Modified Diffuse Reflection Incomplete Accommodation (DRIA) model (Walker et al.)

% Note - Output Cd considers the varying cross-sectional area. No need to
% multiply with area for drag force.

% References - Drag Coefficient Model Using the Cercignani–Lampis–Lord
% Gas–Surface Interaction Model, 10.2514/1.A32677

%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% u = relative velocity vector in the body frame
% Ai = satellite panel normals in body frame (matrix)
% Tw = satellite wall temperature  (K)
% Vi = bulk speed (norm of satellite relative velocity)  (m/s)
% A  = Areas of panels (m^2)
% Ar = Reference area (taken as 1) (m^2)
% M_s = molar mass of satellite surface material for each panel (g/mol)
% Alpha = Energy accommodation coefficient (0 to 1, closer to 1 for diffuse reflection)
% f = fraction of surface covered by atomic oxygen
% M = mean molar mass of the atmospheric gases   (g/mol)
% T = Ambient temperature  (K)
% n_oxa = number density of atomic oxygen (m^-3)
% Kl = Langmuir isotherm parameter
% frac_ox = fraction of surface covered by atomic oxygen (user specified constant or time-varying value calculated here) (0 to 1)

%%%%%%%%%%%%% Outputs %%%%%%%%%%%%%%%%%%%%%%
% Cd_total = total drag-coefficient of satellite

% Constants
R   = 8314.46;           % universal gas constant; Units: g⋅m2⋅s−2⋅K−1⋅mol−1
k_b = 1.38064852e-23;    % Boltzmann constant; Units: m2 kg s-2 K-1

% Calculate input parameters to the model
Ri = R/M;
S  = Vi/sqrt(2*Ri*T);                     % molecular speed ratio


% Fraction of surface covered by atomic oxygen
% If non-physical value specified by user, calculate here.
if frac_ox < 0
    P_o = n_oxa*k_b*T;                 % partial pressure of oxygen (n*kb*T)
    f = Kl.*P_o./(1+Kl.*P_o);
else
    f = frac_ox;
end

% Calculate drag-coefficient
Cd_s = 0;
for j = 1:numel(Ai(1,:))
    Y = u(1)*Ai(1,j) + u(2)*Ai(2,j) + u(3)*Ai(3,j);    % dot product of panel normal and relative velocity vector
    m_r = M/M_s(j);                                          % mass fraction of panels
    Alpha_s = 3.6*m_r*Y/(1+m_r)^2;                           % Goodman's formula for flat plate
    Cd_s(j,:) = sentman(Vi, A(j), Ar,Y,Ri,Alpha_s,Tw,S);
end

Cd_s = sum(Cd_s);
for j = 1:numel(Ai(1,:))
    Y = u(1)*Ai(1,j) + u(2)*Ai(2,j) + u(3)*Ai(3,j);
    Cd(j,:) = sentman(Vi, A(j), Ar,Y,Ri,Alpha,Tw,S);
end
Cd_ads = sum(Cd)+ 0;  %%%%%%%%%%%%%%%%%%%%%%
Cd_total = f*Cd_ads + (1-f)*Cd_s;


% Sentman's formula to calculate the drag-coefficient for each plate
    function Cd = sentman(Vi, A, Ar,Y,Ri,Alpha,Tw,S)
        %% Equations (Doornbos)
        G = 1/(2*S*S);
        P = exp(-Y^2*S^2)/S;
        Q = 1+G;
        Z = 1+erf(Y*S);                %%%%%%% need to find the complement of erf function for Fortran
        Vr = Vi*sqrt(0.5*(1+Alpha*(4*Ri*Tw/Vi^2 -1)));
        
        Cd = (P/sqrt(pi) + Y*Q*Z + Y/2*Vr/Vi*(Y*sqrt(pi)*Z + P))*A/Ar;
    end
end