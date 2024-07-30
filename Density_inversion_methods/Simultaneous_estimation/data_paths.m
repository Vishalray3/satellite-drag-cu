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