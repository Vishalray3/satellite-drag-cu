 
function [lat,lon,h]=ECEF2LLA(x,y,z)

% This function converts x,y,z ECEF position  to geodetic Longitude and Latitude
% and ellipsoidal height WRT to WGS-84 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%%
% x,y,z                      ECEF position in m
%%%%%%%%%%%% Outputs %%%%%%%%%%%%%%
%lat, lon, h                 Geodetic Latitude and longitude respectively 
%                             both in degrees and h is ellipsoidal height
%                             in m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Given/known: 
R=6378.137*10^3; % semimajor axis of WGS-84 ellipsoid in m
f=1/298.257223563; % flattenning parameter of ellipsoid 
esq=2*f-f^2; % square of eccentricity of ellipsoid 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


lon=atan2d(y,x); % Longitude in degrees 
rho=sqrt(x^2+y^2); % projection of user position in xy plane (ECEF coordinates)
r=sqrt(x^2+y^2+z^2); % distance/range to the user location

% To solve for latitude, we need to use numerical methods 
latold=asind(z/r); % Initial guess for geodetic latitude ~= geocentric latitude 
tol=100; % default, can be anything greater than our minimum tolerance 1*10^-8
count=0; % a counter to prevent looping forever or too much 

while tol>= 1*10^-8; % numerically adjusting our latitude guess to within specified tolerance  
    C=R/sqrt(1-esq*(sind(latold))^2);
    lat=atan2d(z+C*esq*sind(latold),rho);
    tol=abs(lat-latold); 
    latold=lat; 
    count=count+1;
    if count>10000 
        break
        error(' Could not converge to within specified tolerance ') 
    end
    end 

% Ellipsoidal height h depends on latitude, so now we can estimate it at this step: 
h=rho/cosd(lat)-C; 



end

