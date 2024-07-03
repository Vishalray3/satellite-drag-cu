%% od plots
load('spire_batch_6_gpst', 'y_res', 'sod_sec_sp3')
figure(1)
plot(sod_sec_sp3/3600, vecnorm(y_res(1:3,:),2,1),'x')
hold on
load('spire_batch_6_gpst_offset', 'y_res', 'sod_sec_sp3')
plot(sod_sec_sp3/3600,vecnorm(y_res(1:3,:),2,1),'x')
xlabel('Time (hrs)')
ylabel('Position residual norm (m)')
legend('POD time in sp3 file','Bias added to POD time in sp3 file')
title('Measurement residuals (position)')
set(gca,'FontSize',14)

load('spire_batch_6_gpst', 'y_res', 'sod_sec_sp3')
figure(2)
plot(sod_sec_sp3/3600, vecnorm(y_res(4:6,:),2,1),'x')
hold on
load('spire_batch_6_gpst_offset', 'y_res', 'sod_sec_sp3')
plot(sod_sec_sp3/3600,vecnorm(y_res(4:6,:),2,1),'x')
xlabel('Time (hrs)')
ylabel('Velocity  residual norm (m/s)')
legend('POD time in sp3 file','Bias added to POD time in sp3 file')
title('Measurement residuals (velocity)')
set(gca,'FontSize',14)


%%
figure(2)
plot(data_.dpre_pos')
ylabel('Prerfit residuals (m)')
title('Prefit residuals with ionospheric-free combination')
set(gca,'FontSize',14)
ylim([-50 50])

figure(3)
plot(data_.freq_diff')
ylabel('L1-L2 (m)')
title('Pseudorange difference between L1-L2')
set(gca,'FontSize',14)