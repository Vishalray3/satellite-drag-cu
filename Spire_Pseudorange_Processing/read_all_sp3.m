function [data_] = read_all_sp3(year,month,dayy,sat_ID, dir_rnx, dir_sp3_q, dir_GPSn)
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
qq = ss; b_sol = []; bv_sol = [];
dsp3_p = []; dsp3_v = []; dsp3_b = []; dsp3_bv = [];
sp3_p = []; sp3_v = []; SS_max_mat = []; SS_min_mat = []; SS_mean_mat = [];



for tm = 1 : num_mon % months loop
    dd_tm = dayy(tm,:); dd_tm=dd_tm(~isnan(dd_tm));
    
    for td = 1 : numel(dd_tm) % days loop
        %%%%%  ############### for loop here depending on the inputs (days/months)
        date = datestr([year,month(tm),dd_tm(td),00,00,00],'yyyy-mm-dd');
        
        
        
        
        %%%% correspinding to this day find files corresponding to a specific satellite ID ...
        %%%%%================= Rinex directory =======================
        S_rnx = dir(dir_rnx);N_rnx = {S_rnx.name}; %  Files in Rinex directory
        ind_date = find(strcmp(N_rnx, date)); %find folder with the specifeid date C=append(dir_rnx,N_rnx{ind_rnx});
        Ss_rnx = dir(append(dir_rnx,N_rnx{ind_date},'/*_0',num2str(sat_ID))); Nn_rnx = {Ss_rnx.name};
        
        
        %%%%%================= sp3 directory =======================
        S_sp3_q = dir(dir_sp3_q);N_sp3_q = {S_sp3_q.name}; %  Files in Rinex directory
        Ss_sp3_q = dir(append(dir_sp3_q,N_sp3_q{ind_date},'/*_0',num2str(sat_ID))); Nn_sp3_q= {Ss_rnx.name};
        
        
        %%%% ################## ANOTHER FOR loop ddepending on size of Nn_rnx
        for j = 1 : numel(Nn_rnx)
            %%%%%%% ================== Read sp3 file ====================
            Sss_sp3=dir([dir_sp3_q N_sp3_q{ind_date} '/' Nn_sp3_q{j} '/*.sp3']);
            filename_sp3 = [Sss_sp3.folder '/' Sss_sp3.name];
            %%%%% Read LEO receiver a-priori ECEF position & velocity coordinates:
            [data, cal_time_gps] = read_pos_sp3(filename_sp3);
            % gps week, time of week in seconds (since 0 hr Sun), sat id, pos, vel
            % (ecef m, m/s), bias in musec
            
            
            
            tsp3 = data(2,:); j = 1;     % seconds of week
%             % Add bias to the seconds of week
%             tsp3 = tsp3 + data(7,:)*1e-6;
            %===================================================================
            %%%% Interpolation of sp3 postions and velocity coordinates:
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
        end % end of multiple files in a day loop
    end % end of days loop
end % end of months loop


data_.recef = recef; data_.vecef = vecef; data_.yyyy=yyyy; data_.mon=mon;
data_.dday=dday; data_.hh=hh; data_.mm=mm; data_.ss=ss; data_.dpost_pos=dpost_p;  data_.dpre_pos=dpre_p;
data_.dpost_vel=dpost_v; data_.dsp3_p=dsp3_p; data_.dsp3_v=dsp3_v;
data_.dsp3_b = dsp3_b; data_.dsp3_bv = dsp3_bv;
data_.qq = qq; data_.b_sol = b_sol; data_.bv_sol = bv_sol;
data_.sp3_p = sp3_p; data_.sp3_v = sp3_v;
data_.ss_max = SS_max_mat; data_.ss_min = SS_min_mat; data_.ss_mean = SS_mean_mat;

%  save('spire_leoAtt_085_.mat','recef','vecef','yyyy','mon','dday','hh','mm','ss','dpost_pos',...
%      'dpost_vel','dsp3_p','dsp3_v');
end




