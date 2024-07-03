%% Re-entry analysis

%% Propagate state and covariance
X_init = [];        % populate from Constants

X_true_aug(:,1) = [X_init(:,1)];
X1 = ode4(@propagator_truth,ode_step,X_true_aug,parameters,Tsamp);
%     X_reg = X1(:,1:6)';
%     X_true = interp1(ode_step, X_reg',time_pred_utc,'spline');
%     X_true = X_true';
%     X_true_full = X1';
%         [t,X1] = ode45(@(t,x)propagator_truth(t, x, parameters), time_prop_utc, X_true_aug(:,1), options);
X_true_aug = X1';
%         deltaV_drag = X_true_aug(7, end);
