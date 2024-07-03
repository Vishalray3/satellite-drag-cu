%% plots k
str_case = 'k_bodf_gpm_case';
ind_case_mc = [1,2,3,4,5,6,7,8,9]; %[1,2,3,5,4];
for ii = 1:numel(ind_case_mc)
    load(strcat(str_case,num2str(ii)))
    k_mean(1,ind_case_mc(ii)) = mean(abs(K_error));
    k_std(1,ind_case_mc(ii)) = std(abs(K_error));
    kvar_mean(1,ind_case_mc(ii)) = mean(sqrt(K_var_final));
    kvar_std(1,ind_case_mc(ii)) = std(sqrt(K_var_final));   
end

vec_box = 2:6;
figure(1)
subplot(2,1,1)
errorbar(k_mean(1,:)/1e6*100,k_std(1,:)/1e6*100,'k*','LineWidth',3)
% errorbar(cd_mean(3,:),cd_std(3,:),'gx','LineWidth',3)
% ylabel('Langmuir constant error (%)')
% xlabel('Case')
% set(gca,'XTickLabel', {'Case 1','Case 2','Case 3','Case 4','Case 5','Case 6','Case 7'});
title('Langmuir constant error (%)')
set(gca,'FontSize',18)
grid on
set(gca,'YScale','log')
yticks([1,10,100])

subplot(2,1,2)
errorbar(kvar_mean(1,:)/1e6*100,kvar_std(1,:)/1e6*100,'k*','LineWidth',3)
% errorbar(cd_mean(3,:),cd_std(3,:),'gx','LineWidth',3)
% ylabel('Standard-deviation (%)')
xlabel('Case')
% set(gca,'XTickLabel', {'Case 1','Case 2','Case 3','Case 4','Case 5','Case 6','Case 7'});
title('Langmuir constant standard-deviation (%)')
set(gca,'FontSize',18)
grid on