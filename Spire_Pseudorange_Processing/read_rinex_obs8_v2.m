function [ rinex ] = read_rinex_obs8_v2(fname, PRN_list, max_epochs)
%*******************************************************
% function [ rinex ] = read_rinex_obs(fname, PRN_list)
%
% DESCRIPTION:
%  
%     This function reads a RINEX format GPS data
%     file and returns the data in an array.
%  
% ARGUMENTS:
%  
%     fname (str) - RINEX file
%     PRN_list (opt) - vector with PRNs to process, useful 
%                      for ignoring some PRN data.
%     epoch_count (Opt) - requested number of epochs
%  
% OUTPUT:
%  
%     rinex - rinex data 
%
% CALLED BY:
%
%
% FUNCTIONS INCLUDED BELOW
%
%     check_rinex_line_length
%     read_rinex_header
%     define_cols
%     cal2jd
%     cal2gps
%     jd2gps
%     adjustyear
%      
%
% MODIFICATIONS:    
% 
%     XX-XX-03  :  Jan Weiss - Original
%     03-14-06  :  Jan Weiss - Cleanup
%               :  See SVN log for further modifications
%     10-26-09  :  P. Axelrad to get rid of conversion of phase to m and
%     validity checking for ASEN 5090
%     10-14-11  :  P. Axelrad - put in check for GPS satellites only
%     10-17-11  :  P. Axelrad - corrected time conversion call
%               :  put in check for comment lines and added C2 data type
%     11-2-12   :  P. Axelrad - now correctly handles files with more than 12
%                  satellites tracked and Glonass, and speeds up reading of 
%                  large files.
%     10-12-13  :  K. Larson, event flag is being read in the wrong place
%                  made a fix, but the logic should be redone. I also 
%                  recommend having epochs sent as the 3rd input instead of 
%                  maxlines. So if someone wants 100 epochs, they can get
%                  it.
%     10-10-14 :   Changed maxlines to epoch count, default is 1 day at 1s
%     10-19-14 :   P. Axelrad Updated to deal with more than 10 observables & to
%                  display each epoch read
%     11-30-14 :   P. Axelrad - updated to deal with arbitrary number of
%                  observations
%     09-25-18 :   P. Axelrad - removed dependence on extra time
%     conversions
%     11-01-18 :   P. Axelrad - fixed problem when reading lots of
%     satellites
%     05-15-21 :   Suood Alnaqbi - modified to read rinex 3.02 for SPIRE
%                  data (CSML Internship CU Boulder - Summer 2021) 
% 
% Colorado Center for Astrodynamics Research
% Copyright 2014 University of Colorado, Boulder
%*******************************************************
blocksize = 100000;

% Initialize variables
rinex_data = [];
line_count = 1;
epoch_count = 0;


% Read header

[ fid, sat_ID, observables ] = read_rinex_header(fname);
num_obs = length(observables);

% Reserve array
r_data = zeros(blocksize,5+num_obs);

icount = 1;


% Status
disp([ 'Parsing RINEX file ' fname ]);
if nargin >= 2 & PRN_list
    disp([ 'Returning data for PRNs: ' num2str(PRN_list) ]);
else
    disp('Returning data for all PRNs');
end

if nargin < 3
    max_epochs = 86400;
end
    
% Get the first line of the observations.
current_line = fgetl(fid);
    
% If not at the end of the file, search for the desired information.
while (current_line ~= -1) & (epoch_count < max_epochs)
    
    
    %Check if this line is a comment line.
    if isempty(strfind(current_line,'COMMENT') )
    
    
    % Check event flag
    event_flag = str2num(current_line(30:32));
    % note from Kristine
    %  the previous line is NOT a good way to find an event_flag
    % an event flag can only come when an epoch begins
    % so that it doesn't crash on existing RINEX files, I 
    % am now also requiring that the first character be blank
    event_flag_more_strict = current_line(30);
    
    if event_flag == 0 & strmatch(event_flag_more_strict,' ');
    
    yr = adjustyear(str2num(current_line(5:6)));
     
    % Get the time for this data epoch.
    current_time = [ yr ; str2num(current_line(8:9)) ; ...
            str2num(current_line(11:12)) ; str2num(current_line(14:15)) ; ...
            str2num(current_line(17:18)) ; str2num(current_line(20:29)) ]';
        
    [gpswk, tow0] = cal2gps(current_time(1:3));
    gpssec = tow0 + current_time(4)*3600 + current_time(5)*60 + current_time(6);
    
    tJD = cal2jd(current_time(1),current_time(2),current_time(3));
    tJD = tJD + current_time(4)*(1/24) + current_time(5)*(1/(24*60)) + current_time(6)*(1/(24*60^2));
    
    
%   [~, gpssec, gpswk] = gpsvec2gpstow(current_time); % older version
    epoch_count = epoch_count+1;
    if rem(epoch_count,100) == 0
    disp([ 'Read ', num2str(epoch_count) ' epochs' ]);
    end

        
    % How many SV's are there?
    current_num_sv = str2num(current_line(34:36));
    current_prn=zeros(current_num_sv,1);
    sat_type=repmat('X',current_num_sv,1);
    
    % Figure out which PRN's there are.
    
    num_sat_lines = ceil(current_num_sv/12)-1; % read extra line(s) if too many satellite for one line
    for i = 1:num_sat_lines
        temp = fgetl(fid);
        current_line = [current_line(1:end) temp(33:min(68,end))];
    end
%     for ii=1:current_num_sv 
%         sat_type(ii) = current_line(30+3*ii);
%         current_prn(ii) = str2num(current_line(31+3*ii : 32+3*ii));        
%     end
%     
    
    % Get the data for all SV's in this epoch.
    for ii=1:current_num_sv
                
        % Get the next line.
        current_line = fgetl(fid);
        [ current_line ] = check_rinex_line_length(current_line);
        line_count = line_count + 1;
        
        sat_type(ii) = current_line(1);
        current_prn(ii) = str2num(current_line(2:3));
        
        

        current_obs = [ str2num(current_line(6:20)); str2num(current_line(22:36)); ...
                        str2num(current_line(38:52));str2num(current_line(54:68)); ... 
                        str2num(current_line(70:84)); str2num(current_line(86:100));...
                         str2num(current_line(102:116)); str2num(current_line(118:131))];  
        obs_count = length(current_obs);
        while obs_count < num_obs            
        current_line = fgetl(fid);
        [ current_line ] = check_rinex_line_length(current_line);
        current_obs = [ current_obs; ...
                        str2num(current_line(6:20)); str2num(current_line(22:36)); ...
                        str2num(current_line(38:52));str2num(current_line(54:68)); ... 
                        str2num(current_line(70:84)); str2num(current_line(86:100));...
                         str2num(current_line(102:116)); str2num(current_line(118:131))];
        obs_count = length(current_obs);
        end
        current_obs = current_obs(1:num_obs);
        
       
       
       
        % Glonass add 100 to PRN, Galileo add 200 to PRN, Compass add 300
        % to PRN
        switch sat_type(ii)
            case {'G',' '} 
                
            case 'R'
                current_prn(ii)=current_prn(ii)+100;
            case 'E'
                current_prn(ii)=current_prn(ii)+200;
            case 'C'
                current_prn(ii)=current_prn(ii)+300;
            case 'S'
                current_prn(ii)=current_prn(ii)+400;
            otherwise
                current_prn(ii)=0;
                continue
        end
        
        % Format the data for this PRN as Date/Time, PRN, Observations.
        current_data = [ gpswk, gpssec, current_prn(ii), sat_ID,tJD, current_obs'];
        
        % Keep only data for the specified PRNs
        if nargin >= 2 & PRN_list & isempty(find(PRN_list == current_prn(ii)))
            continue
        end
           
       
        %Append to the master rinex data file.
        %rinex_data = [ rinex_data ; current_data ];
        r_data(icount,:)= current_data;
        if icount == blocksize
            rinex_data = [rinex_data; r_data];
            r_data = zeros(size(r_data));
            icount = 1;
        else
            icount = icount+1;
        end

        
        
       
    end  % for ii=1:current_num_sv
   
    end % for event flag
    end % for comment line
    % Get the next line.
    current_line = fgetl(fid);
    line_count = line_count + 1;
    
   
end  % while current_line ~= -1
rinex_data = [rinex_data; r_data];
i = find(rinex_data(:,2) == 0);
rinex_data(i,:) =[];

size(rinex_data)
rinex.data = rinex_data;
clear rinex_data

% Define columns
rinex = define_cols(rinex, observables);

% Get rid of zero data entries introduced by
% line padding when more than 5 obs are present
% rinexv3.data = rinexv3.data(:,1:3+num_obs);


% Status
% disp([ 'Total lines: ', num2str(line_count) ]);
disp([ 'Total epochs: ', num2str(epoch_count) ]);
disp('Finished.');
disp(' ');
end

%==========================================================================
function [ current_line ] = check_rinex_line_length(current_line)

if length(current_line) < 131
    
    add_spaces = 131 - length(current_line);
    
    for j = 1 : add_spaces
        
        current_line = [ current_line , '0' ];
        
    end
    
end

% Check if there are any blanks in the data and put a zero there.
current_line = strrep(current_line,' ', '0');
end
%==========================================================================
function [ fid, sat_ID, observables ] = read_rinex_header(file_name )

% Initialize vars
observables = {};
rec_xyz = [ NaN NaN NaN ];
sat_ID = NaN;

% Assign a file ID and open the given header file.
fid=fopen(file_name);

% If the file does not exist, scream bloody murder!
if fid == -1
    display('Error!  Header file does not exist.');
else
    
    % Set up a flag for when the header file is done.
    end_of_header=0;
    
    % Get the first line of the file.
    current_line = fgetl(fid);
    
    % If not at the end of the file, search for the desired information.
    while end_of_header ~= 1
        
        % Search for the approximate receiver location line.
        if strfind(current_line,'APPROX POSITION XYZ')
            
            % Read xyz coordinates into a matrix.
            rec_xyz = sscanf(current_line,'%f');
        end
        
        % Search for SAT ID (MAKER NAME): 
        if strfind(current_line,'MARKER NAME')
            
            % Read xyz coordinates into a matrix.
            sat_ID = sscanf(current_line,'%f');
        end
        
        
        % Search for the number/types of observables line.
        if strfind(current_line,'SYS / # / OBS TYPES')
            
            % Read the non-white space characters into a temp variable.
            [num_obs] = sscanf(current_line(2:end),'%d',1);
            [observables_temp,obs_count] = sscanf(current_line(7:60),'%s');  
            while obs_count < num_obs            
                current_line = fgetl(fid);
                [obs_next,obs_count_next]= sscanf(current_line(7:60),'%s');
                observables_temp = [observables_temp obs_next ];
                obs_count = obs_count+obs_count_next;
            end
                    
            % Read the number of observables space and then create
            % a matrix containing the actual observables.
            for ii = 1:num_obs                
                observables{ii} = observables_temp( 3*ii-2 : 3*ii );
            end
          
        end
        
                  
        % Get the next line of the header file.
        current_line = fgetl(fid);
        
        %Check if this line is at the end of the header file.
        if strfind(current_line,'END OF HEADER')
            end_of_header=1;
        end
        
    end
end
end
%==========================================================================
function rinex = define_cols(rinex, observables)

% Defaults
rinex.col.WEEK = 1; 
rinex.col.TOW = 2;
rinex.col.PRN = 3;
rinex.col.SAT = 4;
rinex.col.JD = 5;

col_offset = 5;

for ii=1:length(observables)

	switch observables{ii}
        case {'L1C'}
            rinex.col.L1C = ii + col_offset;
        case {'L2L'}
            rinex.col.L2L = ii + col_offset;
        case {'LA'}
            rinex.col.LA = ii + col_offset;
        case {'L5'}
            rinex.col.L5 = ii + col_offset;
        case {'P1'}
            rinex.col.P1 = ii + col_offset;
        case {'P2'}
            rinex.col.P2 = ii + col_offset;
        case {'C1C'}
            rinex.col.C1C = ii + col_offset;
        case {'C2L'}
            rinex.col.C2L = ii + col_offset;
        case {'C5'}
            rinex.col.C5 = ii + col_offset;
        case {'S1C'}
            rinex.col.S1C = ii + col_offset;
        case {'S2L'}
            rinex.col.S2L = ii + col_offset;
        case {'SA'}
            rinex.col.SA = ii + col_offset;
        case {'D1C'}
            rinex.col.D1C = ii + col_offset;
        case {'D2L'}
            rinex.col.D2L = ii + col_offset;
    end  % switch
    
end  % for ii=1:length(observables)
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

%==========================================================================

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
function year = adjustyear(yy)
% P. Axelrad - function to correct 2 digit year to the right value

yr=1900*(yy<100); % Convert 00-99 to 1900-1999
yr = yr + 100*(yy<30); % Convert 00-30 to 2000-2030

year = yy+yr; 
end
%==========================================================================



