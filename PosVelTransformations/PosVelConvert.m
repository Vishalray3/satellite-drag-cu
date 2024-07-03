function [r_out,v_out,rot_mat_pos,DUT1,xp,yp] = PosVelConvert(r_in,v_in,EpochUTC,ConversionType,NUTModel, LeapSec, NUT1980Info, EOPInfo)
%
% PosVelConvert    - This function converts position and velocity of a satellite 
%                    between various coordinate frames at some fixed epoch.
%
% Syntax:           [r_out,v_out] = PosVelConvert(r_in,v_in,'EpochUTC','ConversionType','NUTModel')  
%
% Inputs:
%   r_in           - Cartesian position vector (1x3 row vector) [km]
%   v_in           - Cartesian velocity vector (1x3 row vector) [km/s]
%   EpochUTC       - Epoch of the satellite state given in UTC JD secs
%   ConversionType - Choose from the following 30 conversion types:
%                    'ECF2J2K' (Other interchangeable terminology for ECF: ECR,ECEF,ITRF)
%                              (Includes polar motion)
%                    'J2K2ECF' (Other interchangeable terminlogy  for J2K: J2000,MJ2000,EMEJ2000,GCRF,ECIMJ2000)
%                    'TDR2J2K' (TDR is also known as EFG - Excludes polar motion)
%                    'J2K2TDR'
%                    'J2K2TODEarthEquator'
%                    'TODEarthEquator2J2K'
%                    'J2K2MODEarthEquator'
%                    'MODEarthEquator2J2K' 
%                    'J2K2TEME' (Also known as TEME TOD - True of Date)
%                    'TEME2J2K'
%                    'MODEarthEquator2TODEarthEquator'
%                    'TODEarthEquator2MODEarthEquator'
%                    'MODEarthEquator2TDR'
%                    'TDR2MODEarthEquator'
%                    'MODEarthEquator2ECF'
%                    'ECF2MODEarthEquator'
%                    'MODEarthEquator2TEME'
%                    'TEME2MODEarthEquator'
%                    'TODEarthEquator2TDR'
%                    'TDR2TODEarthEquator'
%                    'TODEarthEquator2ECF'
%                    'ECF2TODEarthEquator'
%                    'TODEarthEquator2TEME'
%                    'TEME2TODEarthEquator'
%                    'TDR2ECF'
%                    'ECF2TDR'
%                    'TDR2TEME'
%                    'TEME2TDR'
%                    'ECF2TEME'
%                    'TEME2ECF'
%   NUTModel      -  Choose the number of terms in the nutation model. For
%                    OCM states use 4 terms. [Format: '4terms' or '106terms']
%   LeapSec       -  latest leap seconds
%   NUT1980Info   -  NUT1980 coefficients
%   EOPInfo       -  eop data finals.all IAU1976: for equinox-based
%   transformations
%
% Outputs:
%   r_out         - Cartesian position vector (1x3 row vector) [km]
%   v_out         - Cartesian velocity vector (1x3 row vector) [km/s]
%
% Examples/Validation Cases: 
%
% Other m-files required: DeltaUT1.m, JulianDate.m, LeapSeconds.m, Nutation1980.m, PolarMotion.m
% Subfunctions: None
% MAT-files required: EOP.mat, LeapSeconds.mat, Nutation1980.mat
% Global variables: LeapSecInfo, NUT1980Info, EOPInfo
%
% See also: None
%
% July 2015; Last revision: 13-Jul-2015
%
% ----------------- BEGIN CODE -----------------

    %% =========================== CONSTANTS ==============================
    
    % Earth's angular velocity (rotation rate) [rad/s] (s = solar sec)
    We = 7.2921158553e-5;
    
    % Julian Date - January 1, 2000 12:00 TT (Terrestrial Time)
    JDTT2000 = 2451545.0;
    
    %% ======================== TIME CONVERSIONS ==========================
    
    
    
    % Convert UTC epoch to TAI epoch (different by number of leap seconds)
    
    EpochTAI = EpochUTC + LeapSec;
    
    % Convert TAI epoch to Terrestrial Time (TT) epoch (different by 32.184 sec)
    % datenum/datevec function avoided to eliminate round-off error
    EpochTT = EpochTAI + 32.184;
    
    % Terrestrial Time (TT) - Julian Date
    [JDTT]   = EpochTT/86400;
    
    % Julian centuries (TT) from a particular epoch (i.e. J2000)
    T        = (JDTT - JDTT2000) / 36525;
    
    % Coordinated Universal Time (UTC) and Julian Date (UTC)
    [JDUTC]  = EpochUTC/86400;
    
    % UT1-UTC offset (linearly interpolated)
    [DUT1]   = DeltaUT1(JDUTC, EOPInfo);
    
    % Universal Time (UT1) and Julian Date (UT1)
    EpochUT1 = EpochUTC + DUT1;
    [JDUT1]  = EpochUT1/86400;
    
    % Julian centuries (UT1) from a particular epoch (i.e. J2000)
    TUT1     = (JDUT1 - JDTT2000) / 36525;
    
    
    %% ========================== PRECESSION ==============================
    
    % Precession angles IAU 1976 model. The following values are adopted at epoch J2000.
    % Stats OD book - George Born (pg. 519 - Appendix H.2)
    % Astronomical Almanac (pg. 104 - Section 3.211)
    % Values given in angle seconds (3600" = 1 deg)
    
    zeta  = (2306.2181*T + 0.30188*T^2 + 0.017998*T^3)*(1/3600)*(pi/180);  % radians
    theta = (2004.3109*T - 0.42665*T^2 - 0.041833*T^3)*(1/3600)*(pi/180);  % radians
    z     = (2306.2181*T + 1.09468*T^2 + 0.018203*T^3)*(1/3600)*(pi/180);  % radians

    PREC  = [cos(zeta)*cos(theta)*cos(z)-sin(zeta)*sin(z), -sin(zeta)*cos(theta)*cos(z)-cos(zeta)*sin(z), -sin(theta)*cos(z) ;
             cos(zeta)*cos(theta)*sin(z)+sin(zeta)*cos(z), -sin(zeta)*cos(theta)*sin(z)+cos(zeta)*cos(z), -sin(theta)*sin(z) ;
             cos(zeta)*sin(theta)                        , -sin(zeta)*sin(theta)                        ,  cos(theta)       ];
         
    %% ============================= NUTATION =============================
    
    % Mean longitude of the Moon minus mean longitude of the Moon's perigee
    % (i.e. Moon's mean anomaly) - [Radians]
    L  = (485866.733 + 1717915922.633*T + 31.310*T^2 + 0.064*T^3)*(1/3600)*(pi/180);
    
    % Mean longitude of the Sun minus mean longitude of the Sun's perigee
    % (i.e. Sun's mean anomaly) - [Radians]
    Lp = (1287099.804 + 129596581.224*T - 0.577*T^2 - 0.012*T^3)*(1/3600)*(pi/180);
    
    % Mean longitude of the Moon minus mean longitude of the Moon's node
    % (i.e. mean distance of the Moon from the ascending node) - [Radians]
    F  = (335778.877 + 1739527263.137*T - 13.257*T^2 + 0.011*T^3)*(1/3600)*(pi/180);
    
    % Mean longitude of the Moon minus mean longitude of the Sun
    % (i.e. mean elongation of the Moon from the Sun) - [Radians]
    D  = (1072261.307 + 1602961601.328*T - 6.891*T^2 + 0.019*T^3)*(1/3600)*(pi/180);
    
    % Longitude of the mean ascending node of the lunar orbit on the
    % ecliptic measured from the mean equinox of date - [Radians]
    OM = (450160.28 - 6962890.539*T + 7.455*T^2 + 0.008*T^3)*(1/3600)*(pi/180);
   
    % Mean obliquity of the ecliptic - [Radians]
    mEps = (84381.448 - 46.8150*T - 0.00059*T^2 + 0.001813*T^3)*(1/3600)*(pi/180);
    
    % Obtain nutation terms (Either 4-term or 106-term model)
    [Li,Lpi,Fi,Di,OMi,S1i,S2i,C1i,C2i] = Nutation1980(NUTModel, NUT1980Info);
    
    % Auxilary angle - [Radians]
    Ai   = (Li*L + Lpi*Lp + Fi*F + Di*D + OMi*OM);
    
    % Nutation in obliquity - [Radians]
    dEps = sum(((C1i+C2i*T)*(1/3600)).*cos(Ai))*(pi/180);
    
    % True obliquity of the ecliptic - [Radians]
    Eps  = mEps + dEps;
    
    % Nutation in longitude - [Radians]
    dPsi = sum(((S1i+S2i*T)*(1/3600)).*sin(Ai))*(pi/180);
    
    % The equation of the equinoxes,i.e. nutation in right ascension - [Radians]
    dH   = dPsi * cos(Eps);
    
    % Nutatation Matrix
    NUT  = [    cos(dPsi)     ,              -sin(dPsi)*cos(mEps)               ,              -sin(dPsi)*sin(mEps)                ;
            sin(dPsi)*cos(Eps),  cos(dPsi)*cos(Eps)*cos(mEps)+sin(Eps)*sin(mEps),  cos(dPsi)*cos(Eps)*sin(mEps)-sin(Eps)*cos(mEps) ;
            sin(dPsi)*sin(Eps),  cos(dPsi)*sin(Eps)*cos(mEps)-cos(Eps)*sin(mEps),  cos(dPsi)*sin(Eps)*sin(mEps)+cos(Eps)*cos(mEps)];
        
    %% ========================== SIDERIAL TIME ===========================
    
    % Greenwich Mean Sidereal Time (Vallado, pg. 186, Eq. 3-47)
    GMST = mod(67310.54841 + 3164400184.812866*TUT1 + 0.093104*TUT1^2 - 6.2e-6*TUT1^3,86400)*(1/240)*(pi/180);  % rad (1 sec = 1/240 deg)
        
    % Greenwich Apparent Sidereal Time
    % (Valado uses mean obliquity (mEps) in first correction term,
    % others [Suppl. to Ast. Almanac and Montenbruck) use true obliquity (Eps),
    % the difference is negligible]
    GAST  = GMST + dPsi*cos(Eps) + 0.00264*sin(OM)*(1/3600)*(pi/180) + 0.000063*sin(2*OM)*(1/3600)*(pi/180);
    
    % Siderial Time Matrix
    ST    =      [cos(GAST) , sin(GAST),   0 ;
                  -sin(GAST), cos(GAST),   0 ;
                       0    ,      0   ,   1];
            
    % Time detivative of the ST matrix
    STdot = We * [-sin(GAST), cos(GAST),   0 ;
                  -cos(GAST),-sin(GAST),   0 ;
                       0    ,     0    ,   0];        
    
    %% ========================== POLAR MOTION ============================
   
    % Polar motion constants (linearly interpolated)
    [xp,yp] = PolarMotion(JDUTC, EOPInfo);
    
    % Polar Motion Matrix
    PM      = [cos(xp) , sin(xp)*sin(yp), sin(xp)*cos(yp) ;
                  0    ,     cos(yp)    ,    -sin(yp)     ;
               -sin(xp), cos(xp)*sin(yp), cos(xp)*cos(yp)];
         
    %% ==================== TRANSFORMATION EQUATIONS ======================
    
    switch ConversionType
        
        case 'ECF2J2K'
            % ECF --> ECI J2000 (Position)
            r_out = PREC'*NUT'*ST'*PM'*r_in';
    
            % ECF --> ECI J2000 (Velocity)
            v_out = PREC'*NUT'*STdot'*PM'*r_in' + PREC'*NUT'*ST'*PM'*v_in';
            
            % Rotation matrix 
            rot_mat_pos = PREC'*NUT'*ST'*PM';
            
        case 'J2K2ECF'
            % ECI J2000 --> ECF (Position)
            r_out = PM*ST*NUT*PREC*r_in';
                        
            % ECI J2000 --> ECF (Velocity)
            v_out = PM*STdot*NUT*PREC*r_in' + PM*ST*NUT*PREC*v_in';
            
            % Rotation matrix 
            rot_mat_pos = PM*ST*NUT*PREC;
            
        case 'TDR2J2K'
            % TDR --> ECI J2000 (Position)
            r_out = PREC'*NUT'*ST'*r_in';
    
            % TDR --> ECI J2000 (Velocity)
            v_out = PREC'*NUT'*STdot'*r_in' + PREC'*NUT'*ST'*v_in';
            
            % Rotation matrix 
            rot_mat_pos = PREC'*NUT'*ST';
            
        case 'J2K2TDR'
            % ECI J2000 --> TDR (Position)
            r_out = ST*NUT*PREC*r_in';
                        
            % ECI J2000 --> TDR (Velocity)
            v_out = STdot*NUT*PREC*r_in' + ST*NUT*PREC*v_in';
            
            % Rotation matrix 
            rot_mat_pos = ST*NUT*PREC;
            
        case 'J2K2TODEarthEquator'
            % ECI J2000 --> ECI True of Date Earth Equator (Position)
            r_out = NUT*PREC*r_in';
            
            % ECI J2000 --> ECI True of Date Earth Equator (Velocity)
            v_out = NUT*PREC*v_in';
            
            % Rotation matrix 
            rot_mat_pos = NUT*PREC;
            
        case 'TODEarthEquator2J2K'
            % ECI True of Date Earth Equator --> ECI J2000 (Position)
            r_out = PREC'*NUT'*r_in';
            
            % ECI True of Date Earth Equator --> ECI J2000 (Velocity)
            v_out = PREC'*NUT'*v_in';
            
            % Rotation matrix 
            rot_mat_pos = PREC'*NUT';
            
        case 'J2K2MODEarthEquator'
            % ECI J2000 --> ECI Mean of Date Earth Equator (Position)
            r_out = PREC*r_in';
            
            % ECI J2000 --> ECI Mean of Date Earth Equator (Velocity)
            v_out = PREC*v_in';
            
            % Rotation matrix 
            rot_mat_pos = PREC;
            
        case 'MODEarthEquator2J2K'
            % ECI Mean of Date Earth Equator --> ECI J2000 (Position)
            r_out = PREC'*r_in';
            
            % ECI Mean of Date Earth Equator --> ECI J2000 (Velocity)
            v_out = PREC'*v_in';
            
            % Rotation matrix 
            rot_mat_pos = PREC';
            
        case 'J2K2TEME'
            % ECI J2000 --> ECI TEME (Position)
            r_out = [ cos(dH),sin(dH),0;
                     -sin(dH),cos(dH),0;
                         0   ,   0   ,1]*NUT*PREC*r_in';
            
            % ECI J2000 --> ECI TEME (Velocity)
            v_out = [ cos(dH),sin(dH),0;
                     -sin(dH),cos(dH),0;
                         0   ,   0   ,1]*NUT*PREC*v_in';
                     
            % Rotation matrix 
            rot_mat_pos = [ cos(dH),sin(dH),0;
                            -sin(dH),cos(dH),0;
                                0   ,   0   ,1]*NUT*PREC;         
                        
        case 'TEME2J2K'
            % ECI TEME --> ECI J2000 (Position)
            r_out = PREC'*NUT'*[ cos(dH),sin(dH),0;
                                -sin(dH),cos(dH),0;
                                    0   ,   0   ,1]'*r_in';
            
            % ECI TEME --> ECI J2000 (Velocity)
            v_out = PREC'*NUT'*[ cos(dH),sin(dH),0;
                                -sin(dH),cos(dH),0;
                                    0   ,   0   ,1]'*v_in';
                                
            % Rotation matrix 
            rot_mat_pos = PREC'*NUT'*[ cos(dH),sin(dH),0;
                                        -sin(dH),cos(dH),0;
                                            0   ,   0   ,1]';
                                
        case 'MODEarthEquator2TODEarthEquator'
            % ECI Mean of Date Earth Equator --> ECI True of Date Earth Equator (Position)
            r_out = NUT*r_in';
            
            % ECI Mean of Date Earth Equator --> ECI True of Date Earth Equator (Velocity)
            v_out = NUT*v_in';
            
            % Rotation matrix 
            rot_mat_pos = NUT;
            
        case 'TODEarthEquator2MODEarthEquator'
            % ECI True of Date Earth Equator --> ECI Mean of Date Earth Equator (Position)
            r_out = NUT'*r_in';
            
            % ECI True of Date Earth Equator --> ECI Mean of Date Earth Equator (Velocity)
            v_out = NUT'*v_in';
            
            % Rotation matrix 
            rot_mat_pos = NUT';
            
        case 'MODEarthEquator2TDR'
            % ECI Mean of Date Earth Equator --> TDR (Position)
            r_out = ST*NUT*r_in';
            
            % ECI Mean of Date Earth Equator --> TDR (Velocity)
            v_out = STdot*NUT*r_in' + ST*NUT*v_in';
            
            % Rotation matrix 
            rot_mat_pos = ST*NUT;
            
        case 'TDR2MODEarthEquator'
            % TDR --> ECI Mean of Date Earth Equator (Position)
            r_out = NUT'*ST'*r_in';
            
            % TDR --> ECI Mean of Date Earth Equator (Velocity)
            v_out = NUT'*STdot'*r_in' + NUT'*ST'*v_in';
            
            % Rotation matrix 
            rot_mat_pos = NUT'*ST';
            
        case 'MODEarthEquator2ECF'
            % ECI Mean of Date Earth Equator --> ECF (Position)
            r_out = PM*ST*NUT*r_in';
            
            % ECI Mean of Date Earth Equator --> ECF (Velocity)
            v_out = PM*STdot*NUT*r_in' + PM*ST*NUT*v_in';
            
            % Rotation matrix 
            rot_mat_pos = PM*ST*NUT;
            
        case 'ECF2MODEarthEquator'
            % ECF --> ECI Mean of Date Earth Equator (Position)
            r_out = NUT'*ST'*PM'*r_in';
            
            % ECF --> ECI Mean of Date Earth Equator (Velocity)
            v_out = NUT'*STdot'*PM'*r_in' + NUT'*ST'*PM'*v_in';
            
            % Rotation matrix 
            rot_mat_pos = NUT'*ST'*PM';
            
        case 'MODEarthEquator2TEME'
            % ECI Mean of Date Earth Equator --> ECI TEME (Position)
            r_out = [ cos(dH),sin(dH),0;
                     -sin(dH),cos(dH),0;
                         0   ,   0   ,1]*NUT*r_in';
            
            % ECI Mean of Date Earth Equator --> ECI TEME (Velocity)
            v_out = [ cos(dH),sin(dH),0;
                     -sin(dH),cos(dH),0;
                         0   ,   0   ,1]*NUT*v_in';
                     
            % Rotation matrix 
            rot_mat_pos = [ cos(dH),sin(dH),0;
                            -sin(dH),cos(dH),0;
                                0   ,   0   ,1]*NUT;
                        
        case 'TEME2MODEarthEquator'
            % ECI TEME --> ECI Mean of Date Earth Equator (Position)
            r_out = NUT'*[ cos(dH),sin(dH),0;
                          -sin(dH),cos(dH),0;
                              0   ,   0   ,1]'*r_in';
            
            % ECI TEME --> ECI Mean of Date Earth Equator (Velocity)
            v_out = NUT'*[ cos(dH),sin(dH),0;
                          -sin(dH),cos(dH),0;
                              0   ,   0   ,1]'*v_in';
                          
            % Rotation matrix 
            rot_mat_pos = NUT'*[ cos(dH),sin(dH),0;
                                -sin(dH),cos(dH),0;
                                  0   ,   0   ,1]';
                          
        case 'TODEarthEquator2TDR'
            % ECI True of Date Earth Equator --> TDR (Position)
            r_out = ST*r_in';
            
            % ECI True of Date Earth Equator --> TDR (Velocity)
            v_out = STdot*r_in' + ST*v_in';
            
            % Rotation matrix 
            rot_mat_pos = ST;
            
            
        case 'TDR2TODEarthEquator'
            % TDR --> ECI True of Date Earth Equator (Position)
            r_out = ST'*r_in';
            
            % TDR --> ECI True of Date Earth Equator (Velocity)
            v_out = STdot'*r_in' + ST'*v_in';
            
            % Rotation matrix 
            rot_mat_pos = ST';
            
        case 'TODEarthEquator2ECF'
            % ECI True of Date Earth Equator --> ECF (Position)
            r_out = PM*ST*r_in';
            
            % ECI True of Date Earth Equator --> ECF (Velocity)
            v_out = PM*STdot*r_in' + PM*ST*v_in';
            
            % Rotation matrix 
            rot_mat_pos = PM*ST;
            
        case 'ECF2TODEarthEquator'
            % ECF --> ECI True of Date Earth Equator (Position)
            r_out = ST'*PM'*r_in';
            
            % ECF --> ECI True of Date Earth Equator (Velocity)
            v_out = STdot'*PM'*r_in' + ST'*PM'*v_in';
            
            % Rotation matrix 
            rot_mat_pos = ST'*PM';
            
        case 'TODEarthEquator2TEME'
            % ECI True of Date Earth Equator --> ECI TEME (Position)
            r_out = [ cos(dH),sin(dH),0;
                     -sin(dH),cos(dH),0;
                         0   ,   0   ,1]*r_in';
            
            % ECI True of Date Earth Equator --> ECI TEME (Velocity)
            v_out = [ cos(dH),sin(dH),0;
                     -sin(dH),cos(dH),0;
                         0   ,   0   ,1]*v_in';
                     
            % Rotation matrix 
            rot_mat_pos = [ cos(dH),sin(dH),0;
                            -sin(dH),cos(dH),0;
                                0   ,   0   ,1];
                        
        case 'TEME2TODEarthEquator'
            % ECI TEME --> ECI True of Date Earth Equator (Position)
            r_out = [ cos(dH),sin(dH),0;
                     -sin(dH),cos(dH),0;
                         0   ,   0   ,1]'*r_in';
            
            % ECI TEME --> ECI True of Date Earth Equator (Velocity)
            v_out = [ cos(dH),sin(dH),0;
                     -sin(dH),cos(dH),0;
                         0   ,   0   ,1]'*v_in';
                     
             % Rotation matrix 
            rot_mat_pos = [ cos(dH),sin(dH),0;
                            -sin(dH),cos(dH),0;
                                0   ,   0   ,1]';
                     
        case 'TDR2ECF'
            % TDR --> ECF (Position)
            r_out = PM*r_in';
            
            % TDR --> ECF (Velocity)
            v_out = PM*v_in';
            
            % Rotation matrix 
            rot_mat_pos = PM;
            
        case 'ECF2TDR'
            % ECF --> TDR (Position)
            r_out = PM'*r_in';
            
            % ECF --> TDR (Velocity)
            v_out = PM'*v_in';
            
            % Rotation matrix 
            rot_mat_pos = PM';
            
        case 'TDR2TEME'
            % TDR --> ECI TEME (Position)
            r_out = [ cos(dH),sin(dH),0;
                     -sin(dH),cos(dH),0;
                         0   ,   0   ,1]*ST'*r_in';
                     
            % TDR --> ECI TEME (Velocity)
            v_out = [ cos(dH),sin(dH),0;
                     -sin(dH),cos(dH),0;
                         0   ,   0   ,1]*STdot'*r_in' + [ cos(dH),sin(dH),0;
                                                         -sin(dH),cos(dH),0;
                                                             0   ,   0   ,1]*ST'*v_in';
                 
            % Rotation matrix 
            rot_mat_pos = [ cos(dH),sin(dH),0;
                            -sin(dH),cos(dH),0;
                                0   ,   0   ,1]*ST';
            
        case 'TEME2TDR'
            % TEME --> TDR (Position)
            r_out = ST*[ cos(dH),sin(dH),0;
                        -sin(dH),cos(dH),0;
                            0   ,   0   ,1]'*r_in';
                     
            % TEME --> TDR (Velocity)
            v_out = STdot*[ cos(dH),sin(dH),0;
                           -sin(dH),cos(dH),0;
                               0   ,   0   ,1]'*r_in' + ST*[ cos(dH),sin(dH),0;
                                                            -sin(dH),cos(dH),0;
                                                                0   ,   0   ,1]'*v_in';
            % Rotation matrix 
            rot_mat_pos = ST*[ cos(dH),sin(dH),0;
                            -sin(dH),cos(dH),0;
                                0   ,   0   ,1];
                            
        case 'ECF2TEME'
            % ECF --> ECI TEME (Position)
            r_out = [ cos(dH),sin(dH),0;
                     -sin(dH),cos(dH),0;
                         0   ,   0   ,1]*ST'*PM'*r_in';
                     
            % ECF --> ECI TEME (Velocity)
            v_out = [ cos(dH),sin(dH),0;
                     -sin(dH),cos(dH),0;
                         0   ,   0   ,1]*STdot'*PM'*r_in' + [ cos(dH),sin(dH),0;
                                                             -sin(dH),cos(dH),0;
                                                                 0   ,   0   ,1]*ST'*PM'*v_in';
            % Rotation matrix 
            rot_mat_pos = [ cos(dH),sin(dH),0;
                            -sin(dH),cos(dH),0;
                                0   ,   0   ,1]*ST';                                        
        case 'TEME2ECF'
            % TEME --> ECF (Position)
            r_out = PM*ST*[ cos(dH),sin(dH),0;
                           -sin(dH),cos(dH),0;
                               0   ,   0   ,1]'*r_in';
                     
            % TEME --> ECF (Velocity)
            v_out = PM*STdot*[ cos(dH),sin(dH),0;
                              -sin(dH),cos(dH),0;
                                  0   ,   0   ,1]'*r_in' + PM*ST*[ cos(dH),sin(dH),0;
                                                                  -sin(dH),cos(dH),0;
                                                                      0   ,   0   ,1]'*v_in';
            % Rotation matrix 
            rot_mat_pos = PM*ST*[ cos(dH),sin(dH),0;
                                   -sin(dH),cos(dH),0;
                                    0   ,   0   ,1]';
        otherwise
            error([ConversionType ' as ConversionType is not supported...']);
            
    end
    
    % Transpose output
    r_out = r_out';
    v_out = v_out';
    
return

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
%---------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
%--------------------------------------------------
% D. Plakalovic  | 07-13-2015 |  Re-coded this function from the original 
%                                version (Feb 2013). More efficient, more
%                                functionality.
%             

