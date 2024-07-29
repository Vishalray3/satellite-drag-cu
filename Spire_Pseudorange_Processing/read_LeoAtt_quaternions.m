function [data, cal_time_att, data_eul] = read_LeoAtt_quaternions(fname)
%*******************************************************
% function [ rinex ] = read_rinex_obs(fname, PRN_list)
%
% DESCRIPTION:
%
%     This function reads the quaternioins from LeoAtt
%     file.
%
% ARGUMENTS:
%
%     fname (str) -   file name
%
% OUTPUT:
%
%     data - [GPS week, time into GPS week, satID, x, y, z]
%==========================================================

% Read all lines from the file
data_file = readlines(fname);
num_lines = numel(data_file);

% Preallocate data storage with an estimated size
estimated_lines = 10000; % Change as per expected file size
data = nan(6, estimated_lines); % [gpswk, gpssec, qx, qy, qz, qw]
cal_time_att = nan(estimated_lines, 6); % [year, month, day, hour, minute, second]
data_eul = nan(5, estimated_lines); % [gpswk, gpssec, yaw_deg, pitch_deg, roll_deg]

j = 1;
k = 1;

for i = 1:num_lines
    current_line = data_file{i};
    
    if contains(current_line, '%eof')
        break;
    end
    
    if startsWith(current_line, 'tim')
        % Get the time for this data epoch.
        elems = sscanf(current_line, 'tim %f %f %f %f %f %f %f');
        current_time = elems(1:6)';
        [gpswk, tow0] = cal2gps(current_time(1:3));
        gpssec = tow0 + current_time(4) * 3600 + current_time(5) * 60 + current_time(6);
    end
    
    if startsWith(current_line, 'sca')
        elems = sscanf(current_line, 'sca %f %f %f %f');
        data(:, j) = [gpswk; gpssec; elems];
        cal_time_att(j, :) = current_time;
        j = j + 1;
    end
    
    if startsWith(current_line, 'ang')
        elems = sscanf(current_line, 'ang %f %f %f');
        data_eul(:, k) = [gpswk; gpssec; elems];
        k = k + 1;
    end
end

% Trim preallocated arrays to actual size
data = data(:, 1:j-1);
cal_time_att = cal_time_att(1:j-1, :);
data_eul = data_eul(:, 1:k-1);

end
