%% Linear Kalman Filter
function [x_up, P_up, innov, P_pred,res_post] = NLKF_cc(F, H, Q, R, Omg, y, x0, P0,vec_con)
n_st = length(F);
%% Kalman code

x_pred = F*x0;
P_pred = F*P0*F' + Omg*Q*Omg';
y_pred = H*x_pred;

K = P_pred*H'/(H*P_pred*H' + R);
innov = y - y_pred;
% if innov(3) > pi
%     innov(3) = innov(3) - 2*pi;
% elseif innov(3) < -pi
%     innov(3) = innov(3) + 2*pi;
% end
K(vec_con,:) = 0;
x_up = x_pred + K*innov;

P_up = (eye(n_st) - K*H)*P_pred*(eye(n_st) - K*H)' + K*R*K';

y_up = H*x_up;
res_post = y - y_up;

