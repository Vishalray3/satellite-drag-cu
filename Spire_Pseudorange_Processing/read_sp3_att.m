function [data_] = read_sp3_att(year,month,dayy,sat_ID, dir_sp3_q)
%==================================================
%This code estimates & stores the position and velocity of Spire satellite using
%linearized pseudorange least-squares approach. It reads rinex
%files, GPS navigation, Spire sp3 files & quaternions data corresponding to
% selected time and Spire satellite ID.
%
%%%%Inputs:
% year:               Year number in yyyy format (so far 2018 only)
% month:              Months in a column vector
% dayy:               Days of month corresponding to the month column use nan
%                      to skip a day
% sat_ID:              Spire satellite ID
% pos_iter:           Number of iterations for position solution
% vel_iter:           Number of iterations for velocity solution
% flag_w              Flag, flag_w = 1 for weighted solution flag_w=0 for
%                      unweighted soltuon
%
%%%%Outputs:
% data_               A structural array with the following data:
%                     recef: Spire ECEF position coordinates (m)
%                     vecef: Spire ECEF velocity coordinates (m/s)
%                     b_sol: Spire clock solution (m)
%                     bv_sol: Spire clock rate solution (m/s)
%                     qq: Spire quaternioins values
%                     (yyyy,mon,dday,hh,mm,ss):(year,month,day,hour,minute,seconds)
%                     (dsp3_p,dsp3_v): difference between estimated point
%                                      solution and sp3 results for position (p) in m and
%                                      velocity (v) in m/s.
%
%                     (dsp3_b,dsp3,bv): difference of clock (b) in m and clock rate
%                     soluton (bv) in m/s with sp3 results.
%                     (dpost_p,dpost_v): postfit residuals for positios (p) and velocity (v)
%
%
%   Author: Suood Alnaqbi  -  sual2795@colorado.edu
%   06/17/2021
%  Modification: Vishal Ray - vira0155@colorado.edu
%                1. Rotate LOS vector during time of transmission in ECEF
%                   frame
%                2. Subtracted POD estimated time bias from the RINEX time tags
%                3. Remaining bias in prefits is frequency dependent
%                (different for f1 and f2) and is absorbed in the estimated
%                clock bias

%==================================================


%Given WGS84 GPS constants:
global Omega_E c f1 f2


% %%%% inputs (month,days,satellite ID):
% sat_ID = 85; % Spire satellite ID
% year = 2018; % Year of Spire data files
% month = [9;10]; % Selected Months of Spire data
% dayy = [23,24;1,nan]; % Selected days of Spire data

num_mon = numel(month); % number of months


% pre-allocate to store output:
recef = []; vecef = []; dpost_p=[];  dpre_p=[]; dpost_v = [];
yyyy = []; mon = yyyy; dday  = mon; hh = dday; mm = hh; ss = mm;
yyyy_att = []; mon_att = yyyy; dday_att  = mon; hh_att = dday; mm_att = hh; ss_att = mm;
yyyy_att_all = []; mon_att_all = yyyy; dday_att_all  = mon; hh_att_all = dday; mm_att_all = hh; ss_att_all = mm;

qq = ss; b_sol = []; bv_sol = [];
dsp3_p = []; dsp3_v = []; dsp3_b = []; dsp3_bv = [];
sp3_p = []; sp3_v = []; SS_max_mat = []; SS_min_mat = []; SS_mean_mat = [];
eci_p = []; eci_v = [];

qx_att = []; qy_att = []; qz_att = []; qw_att = [];
qx_att_all = []; qy_att_all = []; qz_att_all = []; qw_att_all = [];


for tm = 1 : num_mon % months loop
    dd_tm = dayy(tm,:); dd_tm=dd_tm(~isnan(dd_tm));
    
    for td = 1 : numel(dd_tm) % days loop
        %%%%%  ############### for loop here depending on the inputs (days/months)
        date = datestr([year,month(tm),dd_tm(td),00,00,00],'yyyy-mm-dd');
        
        
        %%%%%================= sp3 directory =======================
        S_sp3_q = dir(dir_sp3_q);N_sp3_q = {S_sp3_q.name}; %  Files in Rinex directory
        
        N_sp3_date = N_sp3_q(contains(N_sp3_q, date));
        N_sp3_date = N_sp3_date(contains(N_sp3_date, sat_ID));
        %%%% ################## ANOTHER FOR loop ddepending on size of Nn_rnx
        for j = 1 : numel(N_sp3_date)
            %%%%%%% ================== Read sp3 file ====================
            foldername = [dir_sp3_q '/' N_sp3_date{j}];
            Sss_sp3=dir(foldername);
            N_sss_sp3 = {Sss_sp3.name};
            filename_sp3 = N_sss_sp3(contains(N_sss_sp3, 'L1B_leoOrb'));
            if ~isempty(filename_sp3)
                for nn_ind = 1:numel(filename_sp3)
                    path_file_sp3 = [foldername '/' char(filename_sp3{nn_ind}) '/' char(filename_sp3{nn_ind}) ]; % [foldername '/' char(filename_sp3) '/' char(filename_sp3)]; %
                    [data, cal_time_gps] = read_pos_sp3(path_file_sp3);
                    if ~isempty(data)
                        tsp3 = data(2,:);               % seconds of week
                        [tsp3,ind_sp3] = unique(tsp3); data=data(:,ind_sp3); cal_time_gps = cal_time_gps(ind_sp3,:);
                        xr = data(4,:);yr = data(5,:);
                        zr = data(6,:);br = data(7,:);
                        xdr = data(8,:);ydr = data(9,:);
                        zdr = data(10,:);
                        sp3_data = [tsp3;xr;yr;zr;br]; vECEF = [xdr;ydr;zdr];
                        % a-priori LEO receiver ECEF position and velocity coordinates :
                        userECEF = sp3_data(2:4,:)*1000;  %m
                        vuserECEF = vECEF*10^-1; % m/s
                        
                        sp3_p = [sp3_p, userECEF]; sp3_v = [sp3_v, vuserECEF];
                        yyyy = [yyyy, cal_time_gps(:,1)']; mon = [mon, cal_time_gps(:,2)']; dday = [dday, cal_time_gps(:,3)'];
                        hh = [hh, cal_time_gps(:,4)']; mm = [mm, cal_time_gps(:,5)']; ss = [ss, cal_time_gps(:,6)'];
                    end
                end
                
            end
            filename_att = N_sss_sp3(contains(N_sss_sp3, 'L1A_leoAtt'));
            if ~isempty(filename_att)
                for mm_ind = 1:numel(filename_att)
                    
                    path_file_att = [foldername '/' char(filename_att{mm_ind}) '/' char(filename_att{mm_ind})]; %[foldername '/' char(filename_att) '/' char(filename_att)];
                    [data_att,cal_time_att] = read_LeoAtt_quaternions(path_file_att);
                    if ~isempty(data_att)
                        tatt = data_att(2,:);
                        [tatt, ind_att] = unique(tatt); data_att = data_att(:,ind_att);    cal_time_att = cal_time_att(ind_att,:);
                        qx = data_att(3,:);   qy = data_att(4,:);   qz = data_att(5,:);  qw = data_att(6,:);
                        yyyy_att = [yyyy_att, cal_time_att(:,1)']; mon_att = [mon_att, cal_time_att(:,2)']; dday_att = [dday_att, cal_time_att(:,3)'];
                        hh_att = [hh_att, cal_time_att(:,4)']; mm_att = [mm_att, cal_time_att(:,5)']; ss_att = [ss_att, cal_time_att(:,6)'];
                        qx_att = [qx_att, qx];    qy_att = [qy_att, qy];  qz_att = [qz_att, qz];  qw_att = [qw_att, qw];
                    end
                end
            end
            filename_att = N_sss_sp3(contains(N_sss_sp3, 'L1A_telAtt'));
            if ~isempty(filename_att)
                for mm_ind = 1:numel(filename_att)
                    path_file_att = [foldername '/' char(filename_att{mm_ind}) '/' char(filename_att{mm_ind})];
                    [data_att_all,data_eci_all,utc_time_all] = read_highcadence_data(path_file_att);
                    if ~isempty(data_att_all)
                        tatt_all = utc_time_all(:,7);   % unix time
                        [~, ind_att] = unique(tatt_all); data_att_all = data_att_all(:,ind_att); data_eci_all = data_eci_all(:,ind_att); utc_time_all = utc_time_all(ind_att,:);
                        yyyy_att_all = [yyyy_att_all, utc_time_all(:,1)']; mon_att_all = [mon_att_all, utc_time_all(:,2)']; dday_att_all = [dday_att_all, utc_time_all(:,3)'];
                        hh_att_all = [hh_att_all, utc_time_all(:,4)']; mm_att_all = [mm_att_all, utc_time_all(:,5)']; ss_att_all = [ss_att_all, utc_time_all(:,6)'];
                        
                        qx_att_all = [qx_att_all, data_att_all(1,:)];  qy_att_all = [qy_att_all, data_att_all(2,:)]; qz_att_all = [qz_att_all, data_att_all(3,:)];
                        qw_att_all = [qw_att_all, data_att_all(4,:)];
                        eci_p = [eci_p, data_eci_all(1:3,:)]; eci_v = [eci_v, data_eci_all(4:6,:)];
                    end
                end
            end
            
            
        end % end of multiple files in a day loop
    end % end of days loop
end % end of months loop

data_.yyyy=yyyy; data_.mon=mon;
data_.dday=dday; data_.hh=hh; data_.mm=mm; data_.ss=ss;
data_.sp3_p = sp3_p; data_.sp3_v = sp3_v;
data_.qx = qx_att;   data_.qy = qy_att;  data_.qz = qz_att;   data_.qw = qw_att;
data_.yyyy_att=yyyy_att; data_.mon_att=mon_att;
data_.dday_att=dday_att; data_.hh_att=hh_att; data_.mm_att=mm_att; data_.ss_att=ss_att;

data_.yyyy_att_all=yyyy_att_all; data_.mon_att_all=mon_att_all;
data_.dday_att_all=dday_att_all; data_.hh_att_all=hh_att_all; data_.mm_att_all=mm_att_all; data_.ss_att_all=ss_att_all;

data_.eci_p = eci_p; data_.eci_v = eci_v;
data_.qx_all = qx_att_all;   data_.qy_all = qy_att_all;  data_.qz_all = qz_att_all;   data_.qw_all = qw_att_all;
end




