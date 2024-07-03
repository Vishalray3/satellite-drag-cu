%% Read the telAtt spire files
function [data_att_all,data_eci_all,utc_time_all] = read_highcadence_data(path_file_att)

fid = fopen(path_file_att,'r');
if fid == -1
    data_att_all = [];
    data_eci_all = [];
    utc_time_all = [];
else
    
    data_ = readmatrix(path_file_att);
    unix_time = data_(:,2);
    data_att_all = data_(:,3:6)';
    data_eci_all = data_(:,10:15)';
    
    t_date = datetime(unix_time, 'ConvertFrom', 'posixtime');
    utc_time_all(:,1) = year(t_date);
    utc_time_all(:,2) = month(t_date);
    utc_time_all(:,3) = day(t_date);
    utc_time_all(:,4) = hour(t_date);
    utc_time_all(:,5) = minute(t_date);
    utc_time_all(:,6) = second(t_date);
    utc_time_all(:,7) = unix_time;
end
end
