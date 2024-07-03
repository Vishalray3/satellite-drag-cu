%% Get the latest eop data (finals.all (1980)) for equinox-based transformations
% Use finals.all (2000) for CIO-Based transformations
clc
clearvars
year = 2017;

addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Spire_Pseudorange_Processing')
URL = ['https://datacenter.iers.org/data/latestVersion/finals.all.iau1980.txt'];
websave('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder/Density_inversion_methods/data/ancillary_data/finals1.data', URL);

eop_all = read_finals_data(year);
EOPInfo(:,1) = eop_all.fmjd - 29999.5;
EOPInfo(:,2) = eop_all.xp;
EOPInfo(:,3) = eop_all.yp;
EOPInfo(:,4) = eop_all.dut1;

save('EOP', 'EOPInfo')
