%% 
clc
clearvars

var = 'PODnoise_snr_solarmin';
load(var)

if strcmp(var, 'EDRnoise_snr')
arclen_mat = Psig_array;
sigma_pos_mat = sigma_pos_mat(1:5);
load('PODnoise_snr', 'Psig_array')
Psig_array(:,6:7,:,:) = [];
Psig_array(:,:,:,3) = [];
end

if strcmp(var, 'EDRnoise_snr_solarmin')
arclen_mat = Psig_array;
load('PODnoise_snr_solarmin', 'Psig_array')
end
% load('PODnoise_snr')
sz = size(Psig_array);
alt_mat = [300:100:800];
noise_mat = sigma_pos_mat;
for ii=1:sz(1)
    alt = alt_mat(ii);
    for jj = 1:sz(2)
        area = area_mat(jj);
        for kk = 1:sz(3)
            sol_act = solar_level(kk);
            ind = sub2ind(sz(1:3), ii, jj, kk);
            for mm = 1:sz(4)
%                 noisin(mm, ind) = noise_mat(mm);
                altin(mm, ind) = alt;
                areain(mm, ind) = area;
                solin(mm, ind) = sol_act;
                rho_rms(mm, ind) = rho_est_rms_array(ii, jj, kk ,mm);
                rho_mean(mm, ind) = rho_est_mean_array(ii, jj, kk ,mm);
                sig(mm, ind) = Psig_array(ii, jj, kk ,mm);
                
                rho_rms_avg(mm, ind) = rho_est_rms_avg_array(ii, jj, kk ,mm);
                rho_mean_avg(mm, ind) = rho_est_mean_avg_array(ii, jj, kk ,mm);
                
                arc_10(mm, ind) = arc10_array(ii, jj, kk ,mm);
                arc_5(mm, ind) = arc5_array(ii, jj, kk ,mm);
                
%                 if sig(mm, ind) == 0
%                     rho_mean(mm, ind) = NaN;
%                     rho_rms(mm, ind) = NaN;
%                     sig(mm, ind) = NaN;
%                     arc_10(mm, ind) = NaN;
%                     arc_20(mm, ind) = NaN;
%                     rho_rms_avg(mm, ind) = NaN;
%                 end
            end
        end
    end
end


%%
[sig_sort, ind] = sort(sig(1,:));
ind_n = ~isnan(sig_sort);
sig_sort = sig_sort(ind_n);
for ii = 1:numel(rho_rms(:,1))
rho_rms_temp = rho_rms(ii,ind);
rho_rms_sort(ii,:) = rho_rms_temp(ind_n);

area_temp = areain(ii,ind);
area_sort(ii,:) = area_temp(ind_n);

alt_temp = altin(ii,ind);
alt_sort(ii,:) = alt_temp(ind_n);

sol_temp = solin(ii,ind);
sol_sort(ii,:) = sol_temp(ind_n);

arc10_temp = arc_10(ii,ind);
arc5_temp = arc_5(ii,ind);
arc10_sort(ii,:) = arc10_temp(ind_n);
arc5_sort(ii,:) = arc5_temp(ind_n);
rhoavg_temp = rho_rms_avg(ii,ind);
rhoavg_sort(ii,:) = rhoavg_temp(ind_n);
end

if strcmp(var, 'PODnoise_snr')
rho_rms_sort(3,:) = [];
noise_mat(3) = [];
arc10_sort(3,:) = [];
arc5_sort(3,:) = [];
end

edges = [1e-9, 5e-9, 1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4];
% edges = [2.5e-9, 5e-9, 7.5e-9, 1e-8, 2.5e-8, 5e-8, 7.5e-8, 1e-7, 2.5e-7, 5e-7, 7.5e-7, 1e-6, 2.5e-6, 5e-6, 7.5e-6, 1e-5, 2.5e-5, 5e-5, 7.5e-5, 1e-4,...
%     2.5e-4 5e-4, 7.5e-4, 1e-3];
y_ind = discretize(sig_sort, edges);


for ii = 1:numel(edges)-1
    rho_rms_smooth(:,ii) = mean(rho_rms_sort(:,y_ind==ii),2);
    arc10_smooth(:,ii) = mean(arc10_sort(:,y_ind==ii),2)/60;
    arc5_smooth(:,ii) = mean(arc5_sort(:,y_ind==ii),2)/60;
end
% rho_rms_smooth(rho_rms_smooth > 100) = 100;

edges_sm = movmean(edges, 3);
x_grid = repmat(edges(2:end), numel(noise_mat),1);
y_grid = repmat(noise_mat', 1, numel(edges(2:end)));
%%
edges1 = (edges(2:end) + edges(1:end-1))/2;
figure()
% pcolor(x_grid, y_grid, log(rho_rms_smooth))
contourf(edges1, noise_mat*100, (rho_rms_smooth), [0:10:100])
colorbar
xlabel('Drag acceleration ($m/s^2$)',  'interpreter', 'latex')
ylabel('POD noise rms (cm)',  'interpreter', 'latex')
c = colorbar;
c.Label.String = ' Error RMS (%)';
set(gca,'YDir','normal')
title('Relative density error (RMS)');
set(gca,'FontSize',16)
colormap('jet')
set(gca,'XScale','log','Ydir','normal') 
% set(gca,'ColorScale','log')
xticks([1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4])
yticks(noise_mat*100)
grid on

ax = gca;

ax.GridAlpha = 0.8;
hold on

% plot(2.5e-8, 5, 'r.', 'LineWidth', 5,'MarkerSize', 20)
% text(3e-8, 7, 'GDC (AMR: 2.5e-3, H: 375 km, F10.7: 95)', 'Color', 'red', 'FontSize', 14, 'FontWeight', 'bold')


plot(2e-7, 50, 'r.', 'LineWidth', 5,'MarkerSize', 20)
text(1.2e-7, 48, 'Starlink 1', 'Color', 'red', 'FontSize', 14, 'FontWeight', 'bold')   % operational, open-book

plot(2e-6, 50, 'r.', 'LineWidth', 5,'MarkerSize', 20)
text(1.2e-6, 48, 'Starlink 2', 'Color', 'red', 'FontSize', 14, 'FontWeight', 'bold')    % operational, shark-fin

plot(7e-5, 50, 'r.', 'LineWidth', 5,'MarkerSize', 20)
text(1.5e-5, 48, 'Starlink 3', 'Color', 'red', 'FontSize', 14, 'FontWeight', 'bold')    % parking, open-book

plot(3.6e-6, 25, 'r.', 'LineWidth', 5,'MarkerSize', 20)
text(1.8e-6, 23, 'Spire', 'Color', 'red', 'FontSize', 14, 'FontWeight', 'bold')
% 
plot(3.5e-7, 1, 'r.', 'LineWidth', 5,'MarkerSize', 20)
text(2.5e-7, 3, 'GRACE', 'Color', 'red', 'FontSize', 14, 'FontWeight', 'bold')

plot(2.1e-5, 1, 'r.', 'LineWidth', 5,'MarkerSize', 20)
text(1.5e-5, 3, 'GOCE', 'Color', 'red', 'FontSize', 14, 'FontWeight', 'bold')

%%
if strcmp(var, 'EDRnoise_snr')
arc5_smooth(2,3) = 1500;
end
figure()
% pcolor(x_grid, y_grid, log(rho_rms_smooth))
if strcmp(var, 'PODnoise_snr')
contourf(edges1, noise_mat*100, (arc5_smooth), [0, 10, 50, 100, 150, 200, 250, 300, 350, 400, 450])
else
contourf(edges1, noise_mat*100, (arc5_smooth), [0, 10, 50, 100, 200, 300, 400, 500, 600, 700])
end
colorbar
xlabel('Drag acceleration ($m/s^2$)',  'interpreter', 'latex')
ylabel('POD noise rms (cm)',  'interpreter', 'latex')
cbh = colorbar ; %Create Colorbar
cbh.Ticks =  [0, 10, 50, 100, 200, 300, 400, 500, 600, 700];
cbh.Label.String = 'Arc-length (mins.)';
set(gca,'YDir','normal')
title('Arc-length required for density error < 5 %');
set(gca,'FontSize',16)
colormap('jet')
set(gca,'XScale','log','Ydir','normal') 
% set(gca,'ColorScale','log')
xticks([1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3])
yticks(noise_mat*100)
grid on
ax = gca;
ax.GridAlpha = 0.8;

hold on

% plot(2.5e-8, 5, 'r.', 'LineWidth', 5,'MarkerSize', 20)
% text(3e-8, 7, 'GDC (AMR: 2.5e-3, H: 375 km, F10.7: 95)', 'Color', 'red', 'FontSize', 14, 'FontWeight', 'bold')

plot(2e-7, 50, 'r.', 'LineWidth', 5,'MarkerSize', 20)
text(1.2e-7, 48, 'Starlink 1', 'Color', 'red', 'FontSize', 14, 'FontWeight', 'bold')   % operational, open-book

plot(2e-6, 50, 'r.', 'LineWidth', 5,'MarkerSize', 20)
text(1.2e-6, 48, 'Starlink 2', 'Color', 'red', 'FontSize', 14, 'FontWeight', 'bold')    % operational, shark-fin

plot(7e-5, 50, 'r.', 'LineWidth', 5,'MarkerSize', 20)
text(1.5e-5, 48, 'Starlink 3', 'Color', 'red', 'FontSize', 14, 'FontWeight', 'bold')    % parking, open-book


plot(3.6e-6, 25, 'r.', 'LineWidth', 5,'MarkerSize', 20)
text(1.8e-6, 23, 'Spire', 'Color', 'red', 'FontSize', 14, 'FontWeight', 'bold')
% 
plot(3.5e-7, 1, 'r.', 'LineWidth', 5,'MarkerSize', 20)
text(2.5e-7, 3, 'GRACE', 'Color', 'red', 'FontSize', 14, 'FontWeight', 'bold')

plot(2.1e-5, 1, 'r.', 'LineWidth', 5,'MarkerSize', 20)
text(1.5e-5, 3, 'GOCE', 'Color', 'red', 'FontSize', 14, 'FontWeight', 'bold')