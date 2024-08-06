function [parent_directory, dir_data, output_dir] = data_paths(linux_os)

if linux_os == 1
    parent_directory = '/home/vira0155';
    dir_data = '/media/faraday/DATA/thermospheric/spire_data/2022/spire_matlab';
    output_dir = fullfile(dir_data, 'results');
else
    parent_directory = '/Users/vishalray/GoogleDrive/Vishal/PhD/Simulations/Main_working_folder';
    dir_data = fullfile(parent_directory, 'satellite-drag-cu/Density_inversion_methods/Simultaneous_estimation');
    output_dir = dir_data;
end

%% Add path to folders
addpath(fullfile(parent_directory, 'mice/mice/src/mice/'))
addpath(fullfile(parent_directory, 'mice/mice/lib/'))
addpath(fullfile(parent_directory, 'satellite-drag-cu/JB08'))
addpath(fullfile(parent_directory, 'satellite-drag-cu/Density_inversion_methods'))
addpath(fullfile(parent_directory, 'satellite-drag-cu/Density_inversion_methods/data/ancillary_data'))
addpath(fullfile(parent_directory, 'satellite-drag-cu/Density_inversion_methods/data/HASDM_data'))
addpath(fullfile(parent_directory, 'satellite-drag-cu/PosVelTransformations'))
addpath(fullfile(parent_directory, 'satellite-drag-cu/Density_inversion_methods/data/erp_data'))
addpath(fullfile(parent_directory, 'nrlmsis2'))