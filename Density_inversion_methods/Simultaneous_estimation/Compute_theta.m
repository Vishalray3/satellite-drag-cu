%% Compute theta nd rotation matrices
clc
clearvars
load('spire83_propagation_2018_11_7'); %,'jdutc_sec', 'X_true_aug','eop','q_sp3','sod_uni','omega_e','theta_sp3log','sod_sec_sp3log', 'eqeterms')
% Julian Centuries past 1-Jan-2000 12:00 Terrestrial Time
ttt = (jdutc_sec + eop.dat + 32.184 - 86400*2451545.0)/86400/36525; % sec
% Interpolate Delta UT1 to TLE time stamps
dut1 = interp1(86400*eop.fmjd,eop.dut1,jdutc_sec-86400*2400000.5,'linear'); % still in sec
% Julian day in UT1
jdut1 = (jdutc_sec+dut1)/86400;
% Interpolate Earth Orientation parameters (EOP) to TLE time stamps (in UTC)
lod = 1e-3*interp1(86400*eop.fmjd,eop.rlod,jdutc_sec-86400*2400000.5,'linear'); % now in sec
xp = (pi/180/3600)*interp1(86400*eop.fmjd,eop.xp,jdutc_sec-86400*2400000.5,'linear'); % now in radians (from arc-sec)
yp = (pi/180/3600)*interp1(86400*eop.fmjd,eop.yp,jdutc_sec-86400*2400000.5,'linear'); % now in radians (from arc-sec)
dpsi = (pi/180/3600/1000)*interp1(86400*eop.fmjd,eop.dpsi,jdutc_sec-86400*2400000.5,'linear'); % now in radians (from milli-arc-sec)
deps = (pi/180/3600/1000)*interp1(86400*eop.fmjd,eop.deps,jdutc_sec-86400*2400000.5,'linear'); % now in radians (from milli-arc-sec)
% vcorot is for transforming ECEF atmosphere co-rotation v = [0;0;0] and wind = [un;vn;wn] to ECI

for ii = 1:numel(sod_uni)
    [recef,~,~, T_ECI2ECEF] = eci2ecef( X_true_aug(1:3,ii), X_true_aug(4:6,ii), [0;0;0], ttt(ii), jdut1(ii), lod(ii), ...
            xp(ii), yp(ii), eqeterms, dpsi(ii), deps(ii) );
        % corotwind_eci is the velocity of the co-rotating atmosphere (and any wind) relative to the ECI frame
        [~,corotwind_eci,~] = ecef2eci( recef, [0;0;0], ttt(ii), jdut1(ii), lod(ii), ...
            xp(ii), yp(ii), eqeterms, dpsi(ii), deps(ii) );
        
    rot_ECItoVVLH(3,:) = -X_true_aug(1:3,ii)'/norm(X_true_aug(1:3,ii));
    yhat               = cross(X_true_aug(4:6,ii)/norm(X_true_aug(4:6,ii)),X_true_aug(1:3,ii)/norm(X_true_aug(1:3,ii)));
    rot_ECItoVVLH(2,:) = yhat'/norm(yhat);
    rot_ECItoVVLH(1,:) = cross(rot_ECItoVVLH(2,:),rot_ECItoVVLH(3,:));
    rot_ECItoVVLH(1,:) = rot_ECItoVVLH(1,:)/norm(rot_ECItoVVLH(1,:));
%     [~,corotwind_eci,~] = ecef2eci  ( X_true_aug(1:3,ii), [0;0;0],ttt(ii),jdut1(ii),lod(ii),xp(ii),yp(ii),eqeterms,dpsi(ii), deps(ii)  );
    vrel_eci = X_true_aug(4:6,ii) - corotwind_eci;
    vrel_vvlh(:,ii) = rot_ECItoVVLH*vrel_eci;
    
    
    rot_VVLHtoSBF = quat2rotmat(q_sp3(:,ii));
    
    vrel_sbf(:,ii) = rot_VVLHtoSBF*vrel_vvlh(:,ii);
    
    theta_data(ii) = atan2(vrel_sbf(2,ii),vrel_sbf(1,ii));         % angle in x-y plane
    phi_data(ii) = atan2(vrel_sbf(3,ii),sqrt(vrel_sbf(1,ii)^2+vrel_sbf(2,ii)^2));  % angle with z-axis
    
    rot_SBF2ECI_sp3(:,:,ii) = (rot_ECItoVVLH*rot_VVLHtoSBF)';
    rot_ECI2ECEF(:,:,ii) = T_ECI2ECEF;
end
