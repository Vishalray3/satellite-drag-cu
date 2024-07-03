function [bsv,a1,GD] = svclock(ephem_all,t_input, prn)
%==========================================================================
%==========================================================================
% [bsv] = svclock(eph, [week tow], prn)
%
% Calculates the satellite clock corrections from broadcast ephemeris 
%
% Modified by Suood Alnaqbi (ASEN 5090 HW 5) 10/13/2020 to add clock
% corrections calculations and remove previous no longer needed
% functionality.
%
% Author: Ben K. Bradley
% Date: 07/19/2009
%
%
% INPUT:               Description                                  Units
%
%  ephem_all    - matrix of gps satellite orbit parameters           (nx25)
%  
%                  col1: prn, PRN number of satellite
%                  col2: M0, mean anomaly at reference time, rad
%                  col3: delta_n, mean motion difference from computed value, rad/s
%                  col4: ecc, eccentricity of orbit
%                  col5: sqrt_a, square root of semi-major axis, m^0.5
%                  col6: Loa, longitude of ascending node of orbit plane at weekly epoch, rad
%                  col7: incl, inclination angle at reference time, rad
%                  col8: perigee, argument of perigee, rad
%                  col9: ra_rate, rate of change of right ascension, rad/s
%                 col10: i_rate, rate of change of inclination angle, rad/s
%                 col11: Cuc, amplitude of the cosine harmonic correction term to the argument of latitude
%                 col12: Cus, amplitude of the sine harmonic correction term to the argument of latitude
%                 col13: Crc, amplitude of the cosine harmonic correction term to the orbit radius
%                 col14: Crs, amplitude of the sine harmonic correction term to the orbit radius
%                 col15: Cic, amplitude of the cosine harmonic correction term to the angle of inclination
%                 col16: Cis, amplitude of the cosine harmonic correction term to the angle of inclination
%                 col17: Toe, reference time ephemeris (seconds into GPS week)
%                 col18: IODE, issue of data (ephemeris) 
%                 col19: GPS_week, GPS Week Number (to go with Toe)
%                 col20: Toc, time of clock
%                 col21: Af0, satellite clock bias (sec)
%                 col22: Af1, satellite clock drift (sec/sec)
%                 col23: Af2, satellite clock drift rate (sec/sec/sec)
%                 col24: Timing Group Delay (TGD), seconds
%                 col25: health, satellite health (0=good and usable)
%
%  t_input      - GPS times to calculate values at                 [WN TOW] (nx2)
%  prn          - PRN to compute values for (one satellite only)                       
%
%
%
% OUTPUT:       
%    
%  bsv       - satellite clock correction                        (nx1)
%                                     
%
%
% Coupling:
%
%   mean2eccentric.m
%
% References:
% 
%   [1] Interface Control Document: IS-GPS-200D
%         < http://www.navcen.uscg.gov/gps/geninfo/IS-GPS-200D.pdf >
%
%   [2] Zhang, J., et.all. "GPS Satellite Velocity and Acceleration
%         Determination using the Broadcast Ephemeris". The Journal of
%         Navigation. (2006), 59, 293-305.
%            < http://journals.cambridge.org/action/displayAbstract;jsess ...
%                ionid=C6B8C16A69DD7C910989C661BAB15E07.tomcat1?fromPage=online&aid=425362 >
%
%   [3] skyplot.cpp by the National Geodetic Survey
%          < http://www.ngs.noaa.gov/gps-toolbox/skyplot/skyplot.cpp >
%
%
% Last Updated:
%
%  2015/01/22  B.K. Bradley - the capability to look for updated ephem
%                              entries that occur at odd times within each
%                              2hr window has been commented out in this 
%                              function and added to read_GPSbroadcast.m
%                              instead. This moves the computational
%                              overhead to the reading which only occurs
%                              once.
%
%==========================================================================
%==========================================================================


% NOTE: Numbered equations in the code (e.g., Eq. 21) correspond to 
%  equations in the [2] reference.

%==========================================================================
% Load GPS Accepted WGS-84 Constants 
%==========================================================================
global c 
%==========================================================================
% Initialize Output Variables for Speed 
%==========================================================================
sz         = size(t_input,1);
bsv        = ones(sz,1) * NaN;
GD         = ones(sz,1) * NaN; 
a1         = ones(sz,1) * NaN; 


%==========================================================================
% Pull Out Correct Ephemerides 
%==========================================================================

% Pull out ephemerides for PRN in question
kk  = find(ephem_all(:,1) == prn);  % kk is vector containing row numbers of ephem_all that are for sat.no. 'index' 
sat_ephem = ephem_all(kk,:);        % sat_ephem is matrix of all ephem data for each entry of sat.no. 'index'


% No matching PRN found, returning data will be NaNs
if isempty(kk),return,end 



%==========================================================================
% Start Main Calculation Loop 
%==========================================================================

% Compute elapsed times of each ephemeris epoch wrt first entry, seconds
dt_ephem = (sat_ephem(:,19) - sat_ephem(1,19))*604800 + (sat_ephem(:,17) - sat_ephem(1,17));


% Compute elapsed times of each input time wrt first ephemeris entry, seconds
dt_input = (t_input(:,1) - sat_ephem(1,19))*604800 + (t_input(:,2) - sat_ephem(1,17));


for tt = 1:sz % loop through all input times


    % Pull out nearest ephemeris values                                                                                        
    [mn,jj] = min(abs( dt_input(tt) - dt_ephem ));
        
    
                                                      
    if isempty(jj),continue,end  % no matching ephemeris time found. continue to next input time 


    % Pull out common variables from the ephemeris matrix
    %======================================================================
    %toe = sat_ephem(jj,17);           % time of ephemeris
    dt  = dt_input(tt) - dt_ephem(jj); % seconds difference from epoch
    a0   = sat_ephem(jj,21);           % clock bias a0, in sec
    a1(tt,1)   = sat_ephem(jj,22);           % clock drift a1, in sec/sec


    

    % Satellite Clock Correction
    bsv(tt,1) =  c*(a0+a1(tt,1)*dt); 

    % Satellite Group delay correction
    TGD = sat_ephem(jj,24);
    GD(tt,1) = c*TGD;
    



end % END of t_input loop =================================================
%==========================================================================    
