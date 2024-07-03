%% Measurements
function [y_pred, H] = meas_GPS( t, X_nom, meas,N_st)
%% measurement prediction
TEI = cspice_sxform( 'J2000', 'ITRF93', t );


y_pred = TEI*X_nom(1:6);
y_pred = y_pred(meas,1);
%% Jacobian
H = [TEI, zeros(6,N_st-6)];
H = H(meas,:);
end
