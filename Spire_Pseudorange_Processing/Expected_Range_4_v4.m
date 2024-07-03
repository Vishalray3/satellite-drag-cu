function [R,EL,Tt_,count_,satECEF_Tt,vsatECEF_Tt,AZ] = Expected_Range_4_v4(ruserECEF,satECEF_,...
                             prn,gps_ephem,wkn,tow,rcc,bsv,relsv)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function corrects for signal travel time and cacluclates the
% satellite ECEF position and velocity @ time of transmission. 
% - To be improved and organized - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Givens
% Omega_E = 7.2921151467e-5;  % WGS-84 value, rad/s 
% mu =  3.986005e14;     % WGS-84 value, m^3/s^2
% % h = 500*10^3;  % in m
% % RE = 6378.137*10^3;
% % Omega_E = sqrt(mu/(h+RE)^3);

c       = 2.99792458e8;     % GPS accepted speed of light, m/s
tol     = 10^-6;            % Specified tolerance 

% Compute initial uncorrected geometric range for a set of satellites : 
[~,~,RANGE]=compute_azelrange_v2(ruserECEF,satECEF_);
% RANGE = RR; 
% Pre-allocate range,elevation and satellite ECEF position at time of 
% reception : 
R=zeros(1,numel(prn)); EL=R; AZ=R; 
satECEF_Tt=zeros(3,numel(prn));


%ECI to ECEF Rotation: 
%  t0 = tow(1); 
%  theta0 = 0;
%   theta_G = @(t) theta0 + Omega_E*(t - t0);
%    cos_theta_G = @(t) cos(theta_G(t));
%    sin_theta_G = @(t) sin(theta_G(t));
%    A_ECEF_ECI = @(t) [cos_theta_G(t), sin_theta_G(t), 0;...
%                          -sin_theta_G(t), cos_theta_G(t), 0;...
%                             0,  0, 1];
% A_ECI_ECEF = @(t)  [cos_theta_G(t), -sin_theta_G(t), 0;...
%                          sin_theta_G(t), cos_theta_G(t), 0;...
%                             0,  0, 1];

                        

for k=1:numel(RANGE) % LOOP over PRNS 
    
% JD = tJD(k); 
% theta_G = JD2GAST(JD)*pi/180;
%   cos_theta_G = cos(theta_G);
%    sin_theta_G = sin(theta_G);
%    A_ECEF_ECI = [cos_theta_G, sin_theta_G, 0;...
%                          -sin_theta_G, cos_theta_G, 0;...
%                             0,  0, 1];
% A_ECI_ECEF = [cos_theta_G, -sin_theta_G, 0;...
%                          sin_theta_G, cos_theta_G, 0;...
%                             0,  0, 1];

PRN=prn(k);         % satellite PRN 
RANGE_old=RANGE(k); % corresponding uncorrected range 
% bsv_ = bsv(k); relsv_ = relsv(k); 
% Tr=tow - rcc/c + bsv(k)/c + relsv(k)/c; % Time of reception 
Tr=tow - rcc/c; % Time of reception 
% bsv_ = bsv(k);
count=0; % number of iterations counter 
error=1; % initilaize tolerance

%  r0_vec = A_ECI_ECEF*ruserECEF(:,k); 
%  v0_vec = A_ECI_ECEF*vuserECEF(:,k);
%   r0_vec = ruserECEF(:,k); 
%  v0_vec = vuserECEF(:,k);
% %  
% r0=norm(r0_vec); v0=norm(v0_vec); % initial magnitude of velocity and
% %position vectors
% hvec=cross(r0_vec,v0_vec); % Angular momentum vector 
% h=norm(hvec);  hvec_u=hvec/h;
% s_energy=v0^2/2-mu/r0; % specific energy
% a=-mu/(2*s_energy); % semi-major axis
% % ecc_vec=1/mu*(cross(v0_vec,hvec)-mu/r0*r0_vec); ecc=norm(ecc_vec); 
% % ecc_vec_u=ecc_vec/ecc; % eccentrcity 
% n(k) = sqrt(mu/a^3); % mean motion 
% 
% ir_hat = r0_vec/norm(r0_vec); 
% ih_hat = cross(r0_vec,v0_vec)/norm(cross(r0_vec,v0_vec));
% itheta_hat = cross(ih_hat,ir_hat);
% 
% TT = [ir_hat';
%     itheta_hat';
%     ih_hat']; % From ECI to orbital frame 
% TT = TT'; % From orbital to ECI

  while error>=tol && count < 50 % iterate to find range at Tr
      
%     Tt = Tr - RANGE_old'./c + bsv(k)/c + relsv(k)/c; % Time of transmission, in sec 
    Tt = Tr - RANGE_old'./c; % Time of transmission, in sec 

      
    % Satellite ECEF position at time of transmission (Tt)
       [health, pos2_Tt,vel2_Tt] = broadcast_eph2pos_vel(gps_ephem,[wkn Tt],PRN); 
       
       
       %Give ECEF position at Tt for PRN, in m 
       pos2_Tt=pos2_Tt'; 
       vel2_Tt = vel2_Tt';
%        pos2_Tt = A_ECI_ECEF*pos2_Tt;
%        
%        phi=n(k).*(Tr-Tt); % Earth rotation angle 
%         C=TT*[cos(phi),sin(phi),0;
%             -sin(phi),cos(phi),0;
%               0,0,1]*TT'; % Rotation matrix 
%           
%          pos2_Tr=C*(pos2_Tt); % ECI position at time of reception (Tr) 
%          pos2_Tr= A_ECEF_ECI*pos2_Tr; % ECEF position at time of reception (Tr) 

         % compute corrected range, using satECEF at Tr
            if health==0
        [Az,El,Range]=compute_azelrange_v2(ruserECEF,pos2_Tt);
            else 
                Range=NaN; El=NaN; 
                pos2_Tt = nan(3,1); vel2_Tt = nan(3,1);
                Tt_ = nan;
            end
            
        error=norm(Range-RANGE_old); % difference 
        RANGE_old=Range; % reset to iterate 
        
        count=count+1; % number of itereations 
  end
  count_(k) = count;
  Tt_(k) = Tt; 
EL(k)=El; % Elevation angle in deg 
AZ(k)=Az; % Azimuth angle in deg 
R(k)=Range; % expected range in m 
satECEF_Tt(:,k)=pos2_Tt;% satellite ECEF position at time of reception
vsatECEF_Tt(:,k) = vel2_Tt;
end

end
