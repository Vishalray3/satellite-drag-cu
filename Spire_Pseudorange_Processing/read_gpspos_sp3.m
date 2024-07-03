function [data] = read_gpspos_sp3(fname)
%*******************************************************
% function [ rinex ] = read_rinex_obs(fname, PRN_list)
%
% DESCRIPTION:
%
%     This function reads the position from sp3 format
%     file and returns the ECEF xyz position coordinates.
%
% ARGUMENTS:
%
%     fname (str) -  sp3 file name
%
% OUTPUT:
%
%     data - [GPS week, time into GPS week, satID, x, y, z]
%==========================================================

fid = fopen(fname,'r');
j = 1;
file_end = [];

while (isempty(file_end) == 1)
    current_line     = fgetl(fid);
    file_end = findstr(current_line,'EOF');
    
    time_line=findstr(current_line(1),'*');
    if (isempty(time_line) == 0)
        % Get the time for this data epoch.
        current_time = [str2num(current_line(3:7)); str2num(current_line(9:10)) ; ...
            str2num(current_line(12:13)) ; str2num(current_line(15:16)) ; ...
            str2num(current_line(18:19)) ; str2num(current_line(21:end))]';
        
        [gpswk, tow0] = cal2gps(current_time(1:3));
        gpssec = tow0 + current_time(4)*3600 + current_time(5)*60 + current_time(6);
    end
    
    pos_line=findstr(current_line(1),'P');
    if (isempty(pos_line) == 0)
        [char, remain] = strtok(current_line);
        [sat_ID, remain] = strtok(remain);
        [x, remain] = strtok(remain);
        [y, remain] = strtok(remain);
        [z,remain]  = strtok(remain);
        [b]         = strtok(remain);
        
        data_ = [gpswk gpssec str2num(sat_ID) str2num(x) str2num(y) str2num(z) str2num(b)]';
    end
    
    vel_line=findstr(current_line(1),'V');
    if (isempty(vel_line) == 0)
        [char, remain] = strtok(current_line);
        [sat_ID, remain] = strtok(remain);
        [xd, remain] = strtok(remain);
        [yd, remain] = strtok(remain);
        [zd,remain]  = strtok(remain);
        [bd,remain]  = strtok(remain);
        
        data(:,j) = [data_;str2num(xd);str2num(yd);str2num(zd); str2num(bd)];
        j = j + 1;
    end
end



end


%==========================================================================
function [ wn,tow ] = cal2gps(ymd )
% function [ wn,tow ] = cal2gps(ymd )
% INPUT
%   [ymd] = array containing the date as year, month, and day
% OUTPUT
%   wn = Full GPS week number (not modulo 1024)
%   tow =  time of week in seconds
% Notes: intended for use in integer dates not guaranteed to preserve
% sub-second values.
%
jd = cal2jd(ymd(:,1),ymd(:,2),ymd(:,3));
[wn,tow]=jd2gps(jd);
end


function [gpsweek,sow,rollover]=jd2gps(jd);
% JD2GPS  Converts Julian date to GPS week number (since
%   1980.01.06) and seconds of week. Non-vectorized version.
%   See also CAL2JD, DOY2JD, GPS2JD, JD2CAL, JD2DOW, JD2DOY,
%   JD2YR, YR2JD.
% Version: 05 May 2010
% Usage:   [gpsweek,sow,rollover]=jd2gps(jd)
% Input:   jd       - Julian date
% Output:  gpsweek  - GPS week number
%          sow      - seconds of week since 0 hr, Sun.
%          rollover - number of GPS week rollovers (modulus 1024)

% Copyright (c) 2011, Michael R. Craymer
% All rights reserved.
% Email: mike@craymer.com

if nargin ~= 1
    warning('Incorrect number of arguments');
    return;
end
if jd < 0
    warning('Julian date must be greater than or equal to zero');
    return;
end

jdgps = cal2jd(1980,1,6);    % beginning of GPS week numbering
nweek = fix((jd-jdgps)/7);
sow = (jd - (jdgps+nweek*7)) * 3600*24;
rollover = fix(nweek/1024);  % rollover every 1024 weeks
%gpsweek = mod(nweek,1024);
gpsweek = nweek;
end
%==========================================================================
function jd=cal2jd(yr,mn,dy)
% CAL2JD  Converts calendar date to Julian date using algorithm
%   from "Practical Ephemeris Calculations" by Oliver Montenbruck
%   (Springer-Verlag, 1989). Uses astronomical year for B.C. dates
%   (2 BC = -1 yr). Non-vectorized version. See also DOY2JD, GPS2JD,
%   JD2CAL, JD2DOW, JD2DOY, JD2GPS, JD2YR, YR2JD.
% Version: 2011-11-13
% Usage:   jd=cal2jd(yr,mn,dy)
% Input:   yr - calendar year (4-digit including century)
%          mn - calendar month
%          dy - calendar day (including factional day)
% Output:  jd - jJulian date

% Copyright (c) 2011, Michael R. Craymer
% All rights reserved.
% Email: mike@craymer.com

if nargin ~= 3
    warning('Incorrect number of input arguments');
    return;
end
if mn < 1 | mn > 12
    warning('Invalid input month');
    return
end
if dy < 1
    if (mn == 2 & dy > 29) | (any(mn == [3 5 9 11]) & dy > 30) | (dy > 31)
        warning('Invalid input day');
        return
    end
end

if mn > 2
    y = yr;
    m = mn;
else
    y = yr - 1;
    m = mn + 12;
end
date1=4.5+31*(10+12*1582);   % Last day of Julian calendar (1582.10.04 Noon)
date2=15.5+31*(10+12*1582);  % First day of Gregorian calendar (1582.10.15 Noon)
date=dy+31*(mn+12*yr);
if date <= date1
    b = -2;
elseif date >= date2
    b = fix(y/400) - fix(y/100);
else
    warning('Dates between October 5 & 15, 1582 do not exist');
    return;
end
if y > 0
    jd = fix(365.25*y) + fix(30.6001*(m+1)) + b + 1720996.5 + dy;
else
    jd = fix(365.25*y-0.75) + fix(30.6001*(m+1)) + b + 1720996.5 + dy;
end
end
