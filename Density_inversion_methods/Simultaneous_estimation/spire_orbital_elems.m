%% Arrange Spire dataset by orbital elements
clearvars
data_path = '/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/main_working_folder_data/spire/2022_02_05';
addpath(data_path)
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods')
file_names = dir(data_path);
file_names = {file_names.name};

spire_datafiles = file_names(contains(file_names, 'spire_satellite'));

mu_e = 0.3986004415E+15; % m3/s2
rad_e = 6378136.3;

for ii=1:numel(spire_datafiles)
    load(spire_datafiles{ii}, 'Xsp3_eci')
    coe(ii,:) = rv2coe_E(Xsp3_eci(1:3,ii),Xsp3_eci(4:6,ii),mu_e);
    spire_id(ii) = extract(spire_datafiles{ii}, 'FM' + digitsPattern(3));
end
altitude = (coe(:,1) - rad_e)*1e-3;

figure(1)
subplot(2,2,1)
edges_a = [450,455,475:5:490, 505:5:520];
h = histogram(altitude, edges_a);
xlabel('Orbital altitude (km)')
set(gca, 'Fontsize', 12)
xticks(edges_a)
xtickangle(60)


edges = get(h,'BinEdges');

% Create text objects at each location
x = (edges(1:end-1)+edges(2:end))/2;
text(x(3),3,'Icesat-2, 496 km','horizontalalignment','center','verticalalignment','bottom')
xline(496, 'LineWidth', 2)

subplot(2,2,2)
edges_e = [0,5e-4, 1e-3, 1.5e-3, 2e-3, 2.5e-3, 3e-3];
h = histogram(coe(:,2), edges_e);
xlabel('Eccentricity')
set(gca, 'Fontsize', 12)
xticks(edges_e)
edges = get(h,'BinEdges');

% Create text objects at each location
x = (edges(1:end-1)+edges(2:end))/2;
text(x(5),6,'Icesat-2, 0.001398','horizontalalignment','center','verticalalignment','bottom')
xline(0.001398, 'LineWidth', 2)

subplot(2,2,3)
edges_i = [50, 55, 85, 90, 95, 100];
h = histogram(coe(:,3), edges_i);
xlabel('Inclination')
set(gca, 'Fontsize', 12)
xticks(edges_i)
edges = get(h,'BinEdges');
% Create text objects at each location
x = (edges(1:end-1)+edges(2:end))/2;
text(80,10,'Icesat-2, 92','horizontalalignment','center','verticalalignment','bottom')
xline(92, 'LineWidth', 2)

subplot(2,2,4)
edges_i = [-75, -70, 0,5, 40, 45, 85, 90, 95, 100];
histogram(coe(:,4), edges_i)
xlabel('RAAN')
set(gca, 'Fontsize', 12)
% Create text objects at each location
x = (edges(1:end-1)+edges(2:end))/2;
text(30,4,'Icesat-2, 0','horizontalalignment','center','verticalalignment','bottom')
xline(0, 'LineWidth', 2)

sgtitle('Orbital elements of Spire satellites (as of Feb 2022)')


coe(:,1) = altitude;
ind = find(coe(:,1) > 470 & coe(:,1) < 525 & coe(:,2) > 1e-4 & coe(:,2) < 3e-3 & coe(:,3) > 85 & coe(:,3) < 100 & ((coe(:,4) > -10 & coe(:,4) < 10) | (coe(:,4) < -160 | coe(:,4) > 160)));
coe_spire = coe(ind,:);
sat_ids = spire_id(ind);



