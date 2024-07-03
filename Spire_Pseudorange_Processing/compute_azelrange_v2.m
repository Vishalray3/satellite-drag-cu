function [AZ, EL, RANGE] = compute_azelrange_v2(userECEF, satECEF)

% This function computes satellite azimuth, elevation and geometric range.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%
%userECEF                     3 x 1 user ECEF position vector 
%satECEF                      3 x N satellites ECEF position vectors
%%%%%%%%%%%% Outputs%%%%%%%%%%%%%%%
% AZ, EL                    Azimuth and elevation of satellites in degrees
%                           1 x N  vector respectively 
%
% RANGE                     Geometric Range of satellites, in m 
%                           1 x N vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

s=size(satECEF); N=s(2); % number of columns of satECEF array / number of satellites
AZ=zeros(1,N); EL=AZ; RANGE=AZ; % empty arrays to store date for each particular sat
satECEFmat=satECEF; % story the satECEF positions in a matrix

for i=1:N 
    satECEF=satECEFmat(:,i);
LOS_ENU =  compute_LOS_ENU(userECEF, satECEF); % line of sight unit vector in ENU coordinates
AZ(i)=atan2d(LOS_ENU(1),LOS_ENU(2)); % Azimuth in degrees 
EL(i)=asind(LOS_ENU(3)); % Elevation in degrees
RANGE(i)=norm(satECEF-userECEF); % Geometric range in m
end
end
