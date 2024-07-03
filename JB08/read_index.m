%% read indices
clc
clearvars
data_mat = readmatrix('SOLFSMY.txt');  
SOLdata = data_mat(:,2:12)';
save('SOLFSMY', 'SOLdata')
% year, doy, jd, F10, F81c, S10, S81c, M10, M81c, Y10, Y81c

data_mat = readmatrix('DTCFILE.txt');
DTCdata = data_mat(:,2:end)';
save('DTCFILE', 'DTCdata')
%%
geomag_mat = readmatrix('SOLRESAP.txt');
save('SOLRESAP', 'geomag_mat')
% Ap index = 4:11