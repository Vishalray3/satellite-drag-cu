function LOS_ENU =  compute_LOS_ENU(userECEF, satECEF)
% This function computes the normalized line of sight vector from user to
% satellite in East-North-Up (ENU) coordinates 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Suood Alnaqbi - Fall 2020
%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%
%userECEF                     3 x 1 user ECEF position vector 
%satECEF                      3 x 1 satellite ECEF position vector
%%%%%%%%%%%% Outputs%%%%%%%%%%%%%%%
% LOS_ENU                     normalized line of sight vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Line of sight vector (still alligned in ECEF) but origin is
% shifted from Earth's origin to the ENU (user location) as the new origin:
LOS_ECEF=satECEF-userECEF; 

% normalizing the line of sight vector
LOS_ECEF_unitvec=LOS_ECEF/norm(LOS_ECEF);  

% compute the latitude, longitude of user 
[lat,lon]=ECEF2LLA(userECEF(1),userECEF(2),userECEF(3)); 

% compute the transformation matrix from ECEF to ENU
C_ECEF2ENU = ECEF2ENU(lat,lon); 

LOS_ENU=C_ECEF2ENU*LOS_ECEF_unitvec; % The unit line of sight vector in ENU 

end