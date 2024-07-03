function [phi, theta, SBF2ECI] = theta_iner(X_true_aug,X_sun,Rsun,Re,omega_e)
r_hat = X_true_aug(1:3)/norm(X_true_aug(1:3));
v_hat = X_true_aug(4:6)/norm(X_true_aug(4:6));
sun_sat = X_sun ;
sun_sat = sun_sat/norm(sun_sat);

z_hat_body = cross(r_hat,v_hat);
z_hat_body = z_hat_body/norm(z_hat_body);

u1_hat = cross(z_hat_body,sun_sat);
u1_hat = u1_hat/norm(u1_hat);
u2_hat = cross(u1_hat,z_hat_body);
u2_hat = u2_hat/norm(u2_hat);
x_hat_body = u1_hat+u2_hat;
x_hat_body = x_hat_body/norm(x_hat_body);
y_hat_body = u1_hat-u2_hat;
y_hat_body = y_hat_body/norm(y_hat_body);


SBF2ECI = [x_hat_body';y_hat_body';z_hat_body']';
V_body = SBF2ECI'*(X_true_aug(4:6)- cross([0;0;omega_e],X_true_aug(1:3)));
phi = atan2(V_body(3), sqrt(V_body(1)^2+ V_body(2)^2));
theta = atan2(V_body(2),V_body(1));




X = X_true_aug(1:3);               % X_state = [vel;pos] w.r.t to primary
X_rel = X - X_sun;                  %  spacecraft w.r.t third body
a = asin(Rsun/norm(X_rel));
b = asin(Re/norm(X));
c = acos(X'*(X_rel)/(norm(X)*norm(X_rel)));
if c <= abs(a-b)
    fs = 0;
elseif c >= a+b
    fs = 1;
elseif c < a+b && c > abs(a-b)
    x = (c^2+a^2-b^2)/(2*c);
    y = sqrt(a^2-x^2);
    A = a^2*acos(x/a) + b^2*acos((c-x)/b) - c*y;
    fs = 1-A/(pi*a^2);
else
    fs = 90;
end




if fs < 1
    theta = fs*theta;
end
end

