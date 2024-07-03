%%
addpath('/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Density inversion methods/Simultaneous estimation/error analysis')
clc
clearvars
Hmat = [200:5:450]; 
for ii_ind = 1:numel(Hmat)
    ii_ind
    Halt_ind = Hmat(ii_ind);
    Hp_ind = Hmat(ii_ind) ;
    Ha_ind = Hmat(ii_ind) + 5;
    DelV =  trajectory_generator( Halt_ind, Hp_ind, Ha_ind);
    DelV_mat(ii_ind) = DelV;
end
save('DeltaV_dragCancellation_G5')

function [deltaV_drag]= trajectory_generator(Halt_ind, Hp_ind, Ha_ind)
    run Main_simulStudy
end