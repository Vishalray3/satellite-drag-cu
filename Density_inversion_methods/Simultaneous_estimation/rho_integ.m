function Xdot = rho_integ(t,X, parameters)

omega = parameters.omega;
zeta = parameters.zeta;
tau_inv = parameters.tau_inv;
N_st = parameters.N_st;
vec_est = parameters.vec_est;
Ba = parameters.Ba;
Ba = Ba(vec_est);
%% state
Xstate_dot = [0*X(1)+X(2); -omega^2*X(1)-2*zeta*omega*X(2);-tau_inv*X(3)];

%% stm
F = [0, 1,0; -omega^2, -2*zeta*omega,0; 0,0,-tau_inv];


F = F(vec_est,vec_est);
N_ac = numel(vec_est);
stm = reshape(X(N_st+1:N_st+N_ac^2),N_ac,N_ac);
stm_dot = F*stm;

%% noise stm
stm_q_dot = stm*(Ba*Ba')*stm';

%% Integrating the smoothing vector b
b = X(N_st+2*N_ac^2+1:end);
% X_state_dot(vec_est) - F*X(vec_est)
bdot = F*b + Xstate_dot(vec_est) - F*X(vec_est);

Xdot = [Xstate_dot;stm_dot(:);stm_q_dot(:);bdot];





