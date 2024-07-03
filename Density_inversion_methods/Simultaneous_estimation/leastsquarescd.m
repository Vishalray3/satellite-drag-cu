%% Experiment to determine constant Cd 
clc
clearvars
N = 14;
for jj = 1:N
    load(strcat('exp',num2str(jj)))
    Cdn(jj) = parameters.Cd_est;
    rho_corrn(jj,:) = Xs_est(7,:) + Xs_est(9,:);
    Yn(jj,:) = Cdn(jj)*(rho_est+rho_corrn(jj,:));
end

Y = mean(Yn,1)';
R = diag(var(Yn,0,1)); 
Cd_iter = mean(Cd);
rho_iter = mean(rho_corrn,1)';
N_t = numel(Y);
del_x = 0;
for kk = 1:5
    X = [Cd_iter;rho_iter] + del_x;
    H = [rho_est'+X(2:end), X(1)*eye(N_t)];
    Y_pred = X(1)*(rho_est'+X(2:end));
    del_x = H'/(H*H')*(Y - Y_pred);
end
    