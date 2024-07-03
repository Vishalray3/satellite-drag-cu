function [data_, time_sec_sp3, time_sec_att, sod_all, pos_error_true, vel_error_true, sod_error_true] = ...
    remove_data_overlap(data_, day_sat_num, step_size)

ind = ismember(data_.dday, day_sat_num);
data_.yyyy = data_.yyyy(ind);
data_.mon = data_.mon(ind);
data_.dday = data_.dday(ind);
data_.hh = data_.hh(ind);
data_.mm = data_.mm(ind);
data_.ss = data_.ss(ind);
data_.sp3_p = data_.sp3_p(:,ind);
data_.sp3_v = data_.sp3_v(:,ind);

if isfield(data_, 'dday_att')
    ind = ismember(data_.dday_att, day_sat_num);
    data_.yyyy_att = data_.yyyy_att(ind);
    data_.mon_att = data_.mon_att(ind);
    data_.dday_att = data_.dday_att(ind);
    data_.hh_att = data_.hh_att(ind);
    data_.mm_att = data_.mm_att(ind);
    data_.ss_att = data_.ss_att(ind);
    data_.qx = data_.qx(ind);
    data_.qy = data_.qy(ind);
    data_.qz = data_.qz(ind);
    data_.qw = data_.qw(ind);
end

%% data overlaps for POD
time_sec_sp3 = 86400*GREGORIANtoJD_vector(data_.yyyy,data_.mon,data_.dday) + 3600*data_.hh + 60*data_.mm + data_.ss;
ind_over = [];
ind_all = 1:numel(time_sec_sp3);
sod_diff = diff(time_sec_sp3);
ind_diff = find(sod_diff<0);
pos_all = data_.sp3_p;
vel_all = data_.sp3_v;
sod_all = time_sec_sp3 - time_sec_sp3(1);

% pos_overlap_all = double.empty(3,0);
% vel_overlap_all = double.empty(3,0);
ind_over_true = [];
ind_common_true = [];
ind_data_comm_all = [];
for nn = 1:numel(ind_diff)
    ind_inc = find(time_sec_sp3 > time_sec_sp3(ind_diff(nn)));
    if isempty(ind_inc)
        ind_inc = numel(time_sec_sp3) +1;
    end
    ind_temp = [ind_diff(nn)+1: ind_inc(1)-1];
    
    ind_over_start(nn) = ind_temp(1);
    over_interval(nn) = ind_temp(end) - ind_temp(1)+1;
    if over_interval(nn) > ind_over_start(nn)
        over_interval(nn) = ind_over_start(nn) - 1;
        ind_temp = [ind_diff(nn)+1: ind_diff(nn)+over_interval(nn)];
    end
    
    % Beginning of data-arcs
    ind_midpoint = ceil(numel(ind_temp)/2);
    
    ind_over = [ind_over, ind_temp(1:ind_midpoint)];
    
    ind_data_common = [ind_over_start(nn) - over_interval(nn): ind_over_start(nn)-1];
    % End of data-arcs
    ind_data_comm_all = [ind_data_comm_all, ind_data_common(ind_midpoint+1:end)];
    %
    time_data = sod_all(ind_temp);
    time_new = sod_all(ind_data_common);
    %     pos_x = interp1(time_data, pos_all(1, ind_temp), time_new);
    %     pos_y = interp1(time_data, pos_all(2, ind_temp), time_new);
    %     pos_z = interp1(time_data, pos_all(3, ind_temp), time_new);
    %
    %     vel_x = interp1(time_data, vel_all(1, ind_temp), time_new);
    %     vel_y = interp1(time_data, vel_all(2, ind_temp), time_new);
    %     vel_z = interp1(time_data, vel_all(3, ind_temp), time_new);
    %
    %     pos_overlap_all = [pos_overlap_all, [pos_x; pos_y; pos_z]];
    %     vel_overlap_all = [vel_overlap_all, [vel_x; vel_y; vel_z]];
    
    if time_data == time_new
        ind_over_true = [ind_over_true, ind_temp];
        ind_common_true = [ind_common_true, ind_data_common];
    end
    
    
end

% pos_error = pos_all(:, ind_data_comm_all) - pos_overlap_all;
% vel_error = vel_all(:, ind_data_comm_all) - vel_overlap_all;
% sod_error = sod_all(ind_data_comm_all);

pos_error_true = pos_all(:, ind_common_true) - pos_all(:, ind_over_true);
vel_error_true = vel_all(:, ind_common_true) - vel_all(:, ind_over_true);
sod_error_true = sod_all(ind_over_true);

ind_remove_overlap = union(ind_over, ind_data_comm_all);
ind_nonover = setdiff(ind_all, ind_remove_overlap);
data_.yyyy = data_.yyyy(ind_nonover);
data_.mon = data_.mon(ind_nonover);
data_.dday = data_.dday(ind_nonover);
data_.hh = data_.hh(ind_nonover);
data_.mm = data_.mm(ind_nonover);
data_.ss = data_.ss(ind_nonover);
data_.sp3_p = data_.sp3_p(:,ind_nonover);
data_.sp3_v = data_.sp3_v(:,ind_nonover);
time_sec_sp3 = time_sec_sp3(ind_nonover);

%% data overlaps for attitude
if isfield(data_, 'dday_att')
    time_sec_att = 86400*GREGORIANtoJD_vector(data_.yyyy_att,data_.mon_att,data_.dday_att) + 3600*data_.hh_att + 60*data_.mm_att + data_.ss_att;
    
    ind_over_att = [];
    ind_all = 1:numel(time_sec_att);
    sod_diff = diff(time_sec_att);
    ind_diff = find(sod_diff<0);
    for nn = 1:numel(ind_diff)
        ind_inc = find(time_sec_att > time_sec_att(ind_diff(nn)));
        if isempty(ind_inc)
            ind_inc = numel(time_sec_att) +1;
        end
        ind_temp = [ind_diff(nn)+1: ind_inc(1)-1];
        
        ind_over_att = [ind_over_att, ind_temp];
        
    end
    
    ind_nonover_att = setdiff(ind_all, ind_over_att);
    
    data_.yyyy_att = data_.yyyy_att(ind_nonover_att);
    data_.mon_att = data_.mon_att(ind_nonover_att);
    data_.dday_att = data_.dday_att(ind_nonover_att);
    data_.hh_att = data_.hh_att(ind_nonover_att);
    data_.mm_att = data_.mm_att(ind_nonover_att);
    data_.ss_att = data_.ss_att(ind_nonover_att);
    data_.qx = data_.qx(ind_nonover_att);
    data_.qy = data_.qy(ind_nonover_att);
    data_.qz = data_.qz(ind_nonover_att);
    data_.qw = data_.qw(ind_nonover_att);
    
    if isfield(data_, 'yaw')
        data_.yaw = data_.yaw(ind_nonover_att);
        data_.roll = data_.roll(ind_nonover_att);
        data_.pitch = data_.pitch(ind_nonover_att);
    end
    time_sec_att = time_sec_att(ind_nonover_att);
    
    % Remove zero data
    
    ind_remove = (data_.qx == 0) & (data_.qy == 0) & (data_.qz == 0) & (data_.qw == 1);
    
    ind_keep = ~ind_remove;
    
    data_.yyyy_att = data_.yyyy_att(ind_keep);
    data_.mon_att = data_.mon_att(ind_keep);
    data_.dday_att = data_.dday_att(ind_keep);
    data_.hh_att = data_.hh_att(ind_keep);
    data_.mm_att = data_.mm_att(ind_keep);
    data_.ss_att = data_.ss_att(ind_keep);
    data_.qx = data_.qx(ind_keep);
    data_.qy = data_.qy(ind_keep);
    data_.qz = data_.qz(ind_keep);
    data_.qw = data_.qw(ind_keep);
    
    if isfield(data_, 'yaw')
        data_.yaw = data_.yaw(ind_keep);
        data_.roll = data_.roll(ind_keep);
        data_.pitch = data_.pitch(ind_keep);
    end
    time_sec_att = time_sec_att(ind_keep);
else
    time_sec_att = []
end