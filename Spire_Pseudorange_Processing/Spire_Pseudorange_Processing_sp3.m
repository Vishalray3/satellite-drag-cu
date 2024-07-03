function [data_] = Spire_Pseudorange_Processing_sp3(year,month,dayy,sat_ID,pos_iter,vel_iter,flag_w, gps_ephem, dir_rnx, dir_sp3_q, dir_GPSn)
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
            %%%%%%% ================== Read rinex file  ====================
            Sss_rnx=dir([dir_rnx N_rnx{ind_date} '/' Nn_rnx{j} '/*.rnx']);
            filename_rnx = [Sss_rnx.folder '/' Sss_rnx.name];
            rinex = read_rinex_obs8_v2(filename_rnx); % read rinex
            %%%%%%% ================== Read sp3 file ====================
            Sss_sp3=dir([dir_sp3_q N_sp3_q{ind_date} '/' Nn_sp3_q{j} '/*.sp3']);
            filename_sp3 = [Sss_sp3.folder '/' Sss_sp3.name];
            %%%%% Read LEO receiver a-priori ECEF position & velocity coordinates:
            [data] = read_pos_sp3(filename_sp3);
            % gps week, time of week in seconds (since 0 hr Sun), sat id, pos, vel
            % (ecef m, m/s), bias in musec
            %%%%%%% ================== Read quaternions file ====================
            Sss_q=dir([dir_sp3_q N_sp3_q{ind_date} '/' Nn_sp3_q{j} '/*.log']);
            filename_q = [Sss_q.folder '/' Sss_q.name];
            [data_q] = read_LeoAtt_quaternions(filename_q);
            
            %
            % SS=dir(fullfile(C,'*'))
            % X = [SS.isdir] & ~ismember({SS.name},{'.','..'})
            %  NN = {SS(X).name};
            %  [~, ib] = ismember({'2018-09-24T23-00-49_2018-09-24T23-57-15_083','2018-09-24T03-01-19_2018-09-24T03-57-18_083'}, NN)
            %  Y = randperm(numel(NN));
            
            tow = rinex.data(:,rinex.col.TOW); %sec into GPS week
            [tow,ind_u] = unique(tow,'first');trnx = tow;  trnxx=tow;
            wkns = rinex.data(:,rinex.col.WEEK); %GPS week number
            wkns = wkns(ind_u);
            toe=tow(1);        % Time of epoch in sec
            tJD = rinex.data(:,rinex.col.JD); tJD = tJD(ind_u);
            tJDoe_ =tJD(1); %rinex epoch in hrs for plot
            SS = rinex.data(:,rinex.col.S1C); SS_max = max(SS);  % signal strength
            
            
            prn_all = cell(1,numel(tow));
            for j=1:numel(tow)
                epoch=tow(j);
                kk  = find(rinex.data(:,rinex.col.TOW) == epoch); %indexes of each specified epoch
                C1C_ = rinex.data(kk,rinex.col.C1C);  C2L_ = rinex.data(kk,rinex.col.C2L);
                C1C_(C1C_==0) = nan; C2L_(C2L_==0)=nan;
                dC = C1C_-C2L_; dC = dC(~isnan(dC));
                prn_v_ = rinex.data(kk,rinex.col.PRN); prn_v_ = prn_v_(~isnan(dC));
                prn_all{j}=prn_v_; % store observed satellite at each epoch
            end
            
           
            
            tsp3 = data(2,:); j = 1;     % seconds of week
            % Add bias to the seconds of week
            tsp3 = tsp3 + data(7,:)*1e-6;
            %===================================================================
            %%%% Interpolation of sp3 postions and velocity coordinates:
            [tsp3,ind_sp3] = unique(tsp3); data=data(:,ind_sp3);
            for i = 1 : numel(trnx)
                for ii = 1 : numel(tsp3)-1
                    if trnx(i) >= tsp3(ii) && trnx(i) <= tsp3(ii+1) && numel(prn_all{i})>4 % if trnx is within tsp3 incement --> Do interpolation for each value
                        
                        trnx_(j) = trnx(i);
                        k(j) = i;
                        
                        %%% increase tsp3 j counter
                        j = j + 1;
                    end
                end
            end
            trnx = trnx_;
             %%%%%%%%%%%%%% correct the rinex time %%%%%%%%%%%%%%
            br = interp1(tsp3 + data(7,:)*1e-6,data(7,:),trnx); % zeros(1, numel(trnx)); %
            trnx = trnx - br*1e-6;
            tJD = tJD - br*1e-6/86400;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            xr = interp1(tsp3,data(4,:),trnx);yr = interp1(tsp3,data(5,:),trnx);
            zr = interp1(tsp3,data(6,:),trnx); br = interp1(tsp3,data(7,:),trnx); % zeros(1, numel(trnx)); %
            xdr = interp1(tsp3,data(8,:),trnx);ydr = interp1(tsp3,data(9,:),trnx);
            zdr = interp1(tsp3,data(10,:),trnx);
            sp3_data = [trnx;xr;yr;zr;br]; vECEF = [xdr;ydr;zdr];
            % a-priori LEO receiver ECEF position and velocity coordinates :
            userECEF = sp3_data(2:4,:)*1000;  %m
            vuserECEF = vECEF*10^-1; % m/s
            %===================================================================
            
            tq = data_q(2,:); jj = 1; jjj = 1;
            %===================================================================
            %%%% Interpolation+extrapolation of quaternions:
            for ii = 1 : numel(trnx)
                for iii = 1 : numel(tq)-1
                    if trnx(ii) >= tq(iii) && trnx(ii) <= tq(iii+1) % if tq is within tsp3 incement --> Do interpolation for each value
                        
                        t_I(jj) = trnx(ii);
                        k_(jj) = ii;
                        jj = jj + 1;
                    end
                end
            end
            t_E = [trnx(trnx < tq(1)),trnx(trnx > tq(end))];
            qx_I = interp1(tq,data_q(3,:),t_I); qx_E = interp1(tq,data_q(3,:),t_E,'linear','extrap');
            qx = [qx_I,qx_E];
            qy_I = interp1(tq,data_q(4,:),t_I); qy_E = interp1(tq,data_q(4,:),t_E,'linear','extrap');
            qy = [qy_I,qy_E];
            qz_I = interp1(tq,data_q(5,:),t_I); qz_E = interp1(tq,data_q(5,:),t_E,'linear','extrap');
            qz = [qz_I,qz_E];
            qw_I = interp1(tq,data_q(6,:),t_I); qw_E = interp1(tq,data_q(6,:),t_E,'linear','extrap');
            qw = [qw_I,qw_E];
            qs = [qx;qy;qz;qw]; % quaternioins
            %===================================================================
            
            
            
            % Epochs :
            epochs = tow(k); num_epochs=numel(epochs);
            JD = tJD(k);
            
            
            DPr_pre = nan(32,num_epochs); DPr_post = DPr_pre;
            DPrr_pre = DPr_pre; DPrr_post = DPr_pre;
            
            
            % Store epochs data in cell array
            % pre-allocate
            prn_ = cell(1,num_epochs); C1_ = prn_; C2_ = C1_; D_=C1_; S_ = D_;
            for j=1:numel(epochs)
                epoch=epochs(j);
                kk  = find(rinex.data(:,rinex.col.TOW) == epoch); %indexes of each specified epoch
                prn_{j}=rinex.data(kk,rinex.col.PRN); % store observed satellite at each epoch
                C1_{j}=rinex.data(kk,rinex.col.C1C); % store C1 psuedoranges at each epoch
                C2_{j}=rinex.data(kk,rinex.col.C2L); % store C2 psuedoranges at each epoch
                D_{j} = rinex.data(kk,rinex.col.D1C); % store Doppler measurement from L1 signal
                S_{j} = rinex.data(kk,rinex.col.S1C);
            end
            
            % pre-allocate:
            prn_c = nan(32,num_epochs);  D_c = prn_c;  El_c = prn_c;
            satECEF_c = nan(3*num_epochs,32); vsatECEF_c = satECEF_c;
            a1_c = D_c'; S_c = prn_c;
            
            yyyy_ = nan(1,num_epochs); mon_ = yyyy_; day_ = mon_; hh_ = day_; mm_=hh_;
            ss_=mm_;
            
                        %%%%%%%%%%%%%%%%%%%%%%%%%% Change epochs to trnx to reflect that bias has been
            % incorporated in tow%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            epochs = trnx;
            %% Start loop over epochs
            for r=1:num_epochs
                prn=prn_{r}(:); % observed PRNs at this epoch
                
                C1=C1_{r}(:); C2=C2_{r}(:);C1(C1==0)=nan; C2(C2==0)=nan;
                D = D_{r}(:); D(D==0)=nan;
                S = S_{r}(:); S(S==0)=nan;
                
                % % Ionosphere-free psuedorange for each observed sat
                [PIF, ~] = ionocorr(C1, f1, C2, f2);
                ind_PIF = find(~isnan(PIF)); PIF = PIF(ind_PIF);
                prn = prn(ind_PIF); D = D(ind_PIF); S = S(ind_PIF);
                % Satellite ECEF positions at this epoch
                tow=epochs(r); % time of week, in this case an epoch
                wkn=wkns(r); % week number no rollovers - for broadcast ephemeris
                %pre-allocate :
                satECEF=zeros(3,numel(prn)); vsatECEF = satECEF;
                bsv=zeros(1,numel(prn)); relsv=bsv; a1 = bsv;
                GD = relsv;
                
                flag=0; b0=0; b=sp3_data(5,r)*10^-6*c;  iter=0;
                userECEF_ = userECEF(:,r); x=[userECEF_;b];
                
                for i=1:numel(prn)
                    PRN=prn(i);
                    
                    
                    
                    %ECEF position of sat prn # in m
                    % tow_i = tow + (bsv(i)+relsv(i))/c;
                    ind_prn= gps_ephem(3,:) == PRN;
                    gps_states_prn = gps_ephem(:,ind_prn);
                    gps_states = interp1(gps_states_prn(2,:), gps_states_prn(4:11,:)', tow, 'spline'); % km, microsec, dm/s, microsec/sec
                    pos = gps_states(1:3)*1e3; % m
                    vel = gps_states(5:7)*0.1; % m/s
                    
                    %store ECEF positions & velocity of each visible prn as column vectors
                    satECEF(:,i)=pos';
                    vsatECEF(:,i) = vel';
                    
                    R =norm(satECEF(:,i)-userECEF_);
                    counter = 0;
                    error =  1;
                    while and(error>1e-9,counter < 10)
                        Tt(i)=tow-R/c;
                        
                    gps_states = interp1(gps_states_prn(2,:), gps_states_prn(4:11,:)', Tt(i), 'spline'); % km, microsec, dm/s, microsec/sec
                    pos = gps_states(1:3)*1e3; % m
                    vel = gps_states(5:7)*0.1; % m/s
                        
                        satECEF(:,i)=pos';
                        vsatECEF(:,i) = vel';
                        
                        phi=Omega_E*(tow-Tt(i));
                        Rot3=[cos(phi) sin(phi) 0
                            -sin(phi) cos(phi) 0
                            0 0 1];
                        satECEF(:,i)=Rot3*satECEF(:,i);
                        range1=norm(satECEF(:,i)-userECEF_);
                        error=abs(range1-R)/R;
                        R=range1;
                        counter=counter+1;
                    end
                    % Transmitter Clock corrections
                    
                    bsv(i) = gps_states(4)*1e-6*c; % gps clock bias in m
                    a1(i) = gps_states(8)*1e-6*c;
                    % Relativistic corrections
                    relsv(i) = 0; %svrelativity(gps_ephem, [wkn Tt(i)], PRN);
                end
                
                satECEF_Tt = satECEF;
                vsatECEF_Tt = vsatECEF;
                %%%%============== position estimation ===================:
                while flag == 0 && iter <pos_iter
                    
                    %%%% Find geometric range:
                    %                     [~,El,RANGE] = compute_azelrange_v2(userECEF_, satECI_Tt);
                    [~,El,RANGE] = compute_azelrange_v2(userECEF_, satECEF_Tt);
                    
                    
                    
                    ind_El = find(El>-5);
                    if numel(ind_El) >=5
                        PIF = PIF(ind_El); RANGE = RANGE(ind_El);
                        bsv = bsv(ind_El); relsv = relsv(ind_El);
                        prn = prn(ind_El); D = D(ind_El); a1 = a1(ind_El);
                        El = El(ind_El); S = S(ind_El);
                        satECEF_Tt = satECEF_Tt(:,ind_El);
                        satECEF = satECEF(:,ind_El);
                        vsatECEF_Tt = vsatECEF_Tt(:,ind_El);
                        GD = GD(ind_El);
                    end
                    
                    prn_c(prn,r) = prn;
                    D_c(prn,r) = D;
                    S_c(prn,r) = S;
                    a1_c(r,prn') = a1;
                    El_c(prn,r) = El;
                    % [RANGE,~,~,~] = Expected_Range_4(userECEF_,satECEF,...
                    % %                              prn,gps_ephem,wkn,tow,vuserECEF_);
                    % ind_El = find(El > -10); El = El(ind_El); EL{r} = El';
                    % PIF = PIF(ind_El); RANGE = RANGE(ind_El); bsv = bsv(ind_El);
                    % relsv = relsv(ind_El); satECEF = satECEF(:,ind_El);
                    % prn = prn(ind_El);
                    
                    
                    % RANGE_{r} = RANGE; bsv_{r} = bsv; relsv_{r} = relsv;
                    
                    %  dr = satECEF-userECEF_;
                    %  vv = vsatECEF - vuserECEF_;
                    % %  dR = dot(repmat(vuserECEF_,1,numel(prn)),dr)/c;
                    %  dR = dot(vv,dr)/c;
                    
                    
                    % Pre-fit residual
                    % dpr_pre=PIF-(RANGE'-bsv' -relsv' + b);
                    dpr_pre = PIF - (RANGE' -bsv' - relsv' + b);
                    
                    % dpr_pre=PIF-(RANGE');
                    
                    % [dpr_pre,ind_dpr]=rmoutliers(dpr_pre);
                    % satECEF(:,ind_dpr) = [];
                    % prn(ind_dpr) = [];
                    
                    
                    %Construct Signal Strength Weighting Matrix:
                    SS_max = max(S);
                    if flag_w == 1
                        W = diag(S/SS_max);
                    elseif flag_w ==0
                        W = eye(numel(prn));
                    end
                    
                    % Construct A matrix
                    % LOS_ECEF = satECI_Tt-userECI_;
                    LOS_ECEF = satECEF_Tt-userECEF_;
                    eLOS     = LOS_ECEF./(sum(LOS_ECEF.^2)).^0.5;
                    eLOS     =eLOS';
                    G = zeros(numel(prn),4);
                    for j = 1 : numel(prn)
                        G(j,:) = [-eLOS(j,:),1];
                    end
                    Gt=transpose(G);
                    
                    %  Matrix H calculation
                    % H=inv(Gt*G);
                    % Hmat(:,:,r)=H;
                    if iter == 0
                        DPr_pre(prn,r) = dpr_pre;
                        dX=(Gt*W*G)\(Gt*W*dpr_pre);
                        
                    end
                    
                    dX=(Gt*W*G)\(Gt*W*dpr_pre);
                    x=x+dX;
                    % userECI_=x(1:3);
                    % userECEF_ = A_ECEF_ECI(tow - b/c)*userECI_;
                    userECEF_ = x(1:3);
                    % x(1:3) = userECEF_;
                    b=x(4);
                    % dX(1:3) = A_ECEF_ECI(tow - b/c)*dX(1:3);
                    dX_mat(:,r) = dX;
                    
                    
                    if sqrt(sum(dX.^2)) <= 1e-09
                        flag=1;
                    end
                    iter=iter+1;
                end
                
                satECEF_c(3*r-2:3*r,prn') = satECEF_Tt;
                vsatECEF_c(3*r-2:3*r,prn') = vsatECEF_Tt;
                
                % store receiver position+clock solution
                X(:,r) = x;
                % Store SNR parameters
                ss_max(r) = max(S);
                ss_min(r) = min(S);
                ss_mean(r) = mean(S);
                %Geometric Range :
                % [~,~,RANGE] = compute_azelrange_v2(userECEF_, satECEF);
                
                %post-fit residuals
                dpr_post=dpr_pre-G*dX;
                % dpr_post = PIF-(RANGE'-bsv' -relsv' + b);
                % for oo = 1 : numel(El)
                %     if El(oo) < 0
                %         dpr_post(oo) = nan;
                %     end
                % end
                DPr_post(prn,r) = dpr_post;
                clear Tt bsv relsv PIF RANGE bsv relsv prn  D a1 El  S satECEF_Tt satECEF vsatECEF_Tt GD 
            end % end of epochs loop
            
            % apriori receiver clock rate:
            coefficients = polyfit(epochs, X(4,:), 1);bd = coefficients(1);
            
            %%%%============== velocity estimation ===================:
            
            for r = 1 : num_epochs
                
                D = D_c(~isnan(prn_c(:,r)),r);
                S = S_c(~isnan(prn_c(:,r)),r);
                a1 = a1_c(r,~isnan(prn_c(:,r)));
                
                
                userECEF_  = X(1:3,r);
                % userECI_   = A_ECI_ECEF(tow)*userECEF_;
                vuserECEF_ = vuserECEF(:,r);
                % vuserECI_ =  A_ECI_ECEF(tow)*(cross([0;0;Omega_E],userECEF_)+vuserECEF(:,r));
                % vuserECI_h(:,r) = vuserECI_;
                satECEF= satECEF_c(3*r-2:3*r,~isnan(prn_c(:,r))');
                vsatECEF = vsatECEF_c(3*r-2:3*r,~isnan(prn_c(:,r))');
                
                bd_ = bd;
                xd = [vuserECEF_;bd_];
                flag2=0; iter2=0;
                while flag2 == 0 && iter2 <vel_iter
                    
                    % Updated LOS unit vector from position solution:
                    LOS_ECEF = satECEF-userECEF_;
                    eLOS     = LOS_ECEF./(sum(LOS_ECEF.^2)).^0.5;
                    eLOS     =eLOS';
                    
                    %Construct the knowns from the Doppler shift equation for range-rate:
                    rrsat = dot(eLOS',vsatECEF);
                    rruser = dot(eLOS',repmat(vuserECEF_,1,numel(a1)));
                    eta = (1+a1')./(1+(rrsat'/c));
                    dprr_pre = -c/f1*D  - eta.*rrsat' + eta.*rruser' + c*a1' - (1+D/f1)*bd_;
                    
                    DPrr_pre(~isnan(prn_c(:,r)),r) = dprr_pre;
                    
                    
                    %Construct Signal Strength Weighting Matrix:
                    if flag_w == 1
                        W = diag(S/SS_max);
                    elseif flag_w ==0
                        W = eye(numel(a1));
                    end
                    
                    
                    
                    % Construct Gv matrix:
                    Gv = zeros(numel(a1),4);
                    for j = 1 : numel(a1)
                        Gv(j,:) = [-eLOS(j,:)*eta(j),(1+D(j)/f1)];
                    end
                    Gvt = transpose(Gv);
                    
                    dXd=(Gvt*W*Gv)\(Gvt*W*dprr_pre);
                    dXd_mat(:,r) = dXd;
                    xd=xd+dXd;
                    vuserECEF_ = xd(1:3);
                    bd_ = xd(4);
                    
                    if sqrt(sum(dXd.^2)) <= 1e-07
                        flag2=1;
                    end
                    iter2=iter2+1;
                end
                
                
                Xd(:,r) = xd;
                
                %post-fit residuals
                dprr_post=dprr_pre-Gv*dXd;
                % rruser = dot(eLOS',repmat(vuserECEF_,1,numel(a1)));
                % dprr_post = -c/f1*D  - eta.*rrsat' + eta.*rruser' + c*a1' - (1+D/f1)*bd_;
                DPrr_post(~isnan(prn_c(:,r)),r) = dprr_post;
                
                [yyyy_(r),mon_(r),day_(r),hh_(r),mm_(r),ss_(r)] = JDtoGREGORIAN_vector(JD(r));
            end % end of epochs loop
            
            %%%% Store values:
            recef = [recef,X(1:3,:)]; vecef = [vecef,Xd(1:3,:)]; b_sol = [b_sol,X(4,:)];
            bv_sol = [bv_sol,Xd(4,:)];
            yyyy = [yyyy,yyyy_]; mon = [mon,mon_]; dday = [dday, day_];
            hh = [hh,hh_]; mm = [mm,mm_]; ss = [ss,ss_];
            qq = [qq,qs]; dpost_p = [dpost_p,DPr_post]; dpre_p = [dpre_p,DPr_pre]; dpost_v = [dpost_v,DPrr_post];
            dsp3_p = [dsp3_p,userECEF-X(1:3,:)]; dsp3_v = [dsp3_v,vuserECEF-Xd(1:3,:)];
            dsp3_b = [dsp3_b,sp3_data(5,:)*10^-6*c-X(4,:)];
            dsp3_bv = [dsp3_bv,sp3_data(5,:)*10^-6*c-X(4,:)];
            sp3_p = [sp3_p, userECEF]; sp3_v = [sp3_v, vuserECEF];
            SS_max_mat = [SS_max_mat, ss_max]; SS_min_mat = [SS_min_mat, ss_min];
            SS_mean_mat = [SS_mean_mat, ss_mean];
            
            clear trnx_ k t_I k_ dX_mat X dXd_mat Xd ss_max ss_min ss_mean
        end % end of multiple files in a day loop
    end % end of days loop
end % end of months loop


data_.recef = recef; data_.vecef = vecef; data_.yyyy=yyyy; data_.mon=mon;
data_.dday=dday; data_.hh=hh; data_.mm=mm; data_.ss=ss; data_.dpost_pos=dpost_p; data_.dpre_pos=dpre_p;
data_.dpost_vel=dpost_v; data_.dsp3_p=dsp3_p; data_.dsp3_v=dsp3_v;
data_.dsp3_b = dsp3_b; data_.dsp3_bv = dsp3_bv;
data_.qq = qq; data_.b_sol = b_sol; data_.bv_sol = bv_sol;
data_.sp3_p = sp3_p; data_.sp3_v = sp3_v;
data_.ss_max = SS_max_mat; data_.ss_min = SS_min_mat; data_.ss_mean = SS_mean_mat;

%  save('spire_leoAtt_085_.mat','recef','vecef','yyyy','mon','dday','hh','mm','ss','dpost_pos',...
%      'dpost_vel','dsp3_p','dsp3_v');
end




