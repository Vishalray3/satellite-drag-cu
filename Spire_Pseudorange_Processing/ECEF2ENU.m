

function C_ECEF2ENU = ECEF2ENU(lat,lon)

%This function calculates the transformation matrix to rotate from ECEF to
%ENU

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%
%lat, lon                     Geodetic Latitude and longitude respectively 
%                             (both in degrees) 
%%%%%%%%%%%% Outputs%%%%%%%%%%%%%%%
% C_ECEF2ENU                     The 3 x 3 transformation matrix from ECEF
%                                to ENU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%Computing the sine and cosine of latitude and longitude to make
%easier next steps 

sinlon=sind(lon); coslon=cosd(lon);
sinlat=sind(lat); coslat=cosd(lat);

%The transformation matrix C from ECEF to ENU :
C_ECEF2ENU=[-sinlon, coslon, 0;
            -sinlat*coslon, -sinlat*sinlon, coslat;
            coslat*coslon, coslat*sinlon, sinlat];

end