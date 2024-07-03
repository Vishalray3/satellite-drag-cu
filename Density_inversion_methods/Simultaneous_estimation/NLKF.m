%% Linear Kalman Filter
function [x_up,sig, P_up, innov, P_pred,res_post] = NLKF(F, H, Q, R, Omg, y, x0, P0)
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
x_up = x_pred + K*innov;
P_up = (eye(n_st) - K*H)*P_pred*(eye(n_st) - K*H)' + K*R*K';
sig = sqrt(diag(P_up));
y_up = H*x_up;
res_post = y - y_up;
% res_post;
