function [data, cal_time_gps] = read_pos_sp3(fname)
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


% Preallocate data storage with an estimated size
estimated_lines = 10000; % Change as per expected file size
data = nan(10, estimated_lines); % Assuming a fixed number of elements per line
cal_time_gps = nan(estimated_lines, 6); % Assuming 6 time elements


% Read the file content
data_file = readlines(fname);
num_lines = numel(data_file);

j = 1;


% Iterate through each line of the file
for i = 1:num_lines
    current_line = data_file{i};
    
    if startsWith(current_line, 'EOF')
        break;
    end
    
    if startsWith(current_line, '*')
        % Get the time for this data epoch.
        elems = sscanf(current_line, '* %f %f %f %f %f %u');
        current_time = elems(1:6)';
        
        [gpswk, tow0] = cal2gps(current_time(1:3));
        gpssec = tow0 + current_time(4) * 3600 + current_time(5) * 60 + current_time(6);
    end
    
    if startsWith(current_line, 'P')
        elems = sscanf(current_line, 'P%u %f %f %f %f');
        sat_ID = elems(1);

        data_ = [gpswk; gpssec; sat_ID; elems(2:5)];
    end
    
    if startsWith(current_line, 'V')
        elems = sscanf(current_line, 'V%u %f %f %f %f');
        
        % Store the data
        data(:, j) = [data_; elems(2:4)];
        cal_time_gps(j, :) = current_time;
        j = j + 1;
    end
end

% Trim preallocated arrays to actual size
data = data(:, 1:j-1);
cal_time_gps = cal_time_gps(1:j-1, :);

end
