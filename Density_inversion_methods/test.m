clearvars
clc
for ind_time = 1:500
    ind_time
    del_T = 10;
Tlast = ind_time*100;
time_prop_utc = 0:del_T:Tlast;
run Main_obsStudy
state_err(:,ind_time) = X_true_aug(1:N_st) - X_nom(1:N_st);
clearvars -except state_err ind_time del_T
end