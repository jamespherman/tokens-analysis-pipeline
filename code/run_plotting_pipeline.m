% run_plotting_pipeline.m
%
% Description:
%   This script serves as the main entry point for generating all final,
%   aggregated figures for the project. It loads aggregated data and then
%   calls a series of plotting functions to generate the figures.
%
% Author:
%   Jules
%
% Date:
%   2025-09-12

% Start timer and provide user feedback
tic;
disp('Starting the plotting pipeline...');

%% Setup Paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);
aggFileName = fullfile(project_root, 'data', 'processed', ...
    'aggregated_analysis_data.mat');

% Load aggregated data
disp('Loading aggregated analysis data...');
data = load(aggFileName);
disp('Data loaded successfully.');

% Define the data types to process
data_types = {'sc', 'snc'};

for i = 1:numel(data_types)
    type = data_types{i};
    data_field = ['aggregated_' type '_data'];
    label = upper(type);

    if isfield(data, data_field)
        current_data = data.(data_field);

        disp(['--- Generating ' label ' Plots ---']);

        disp(['Generating aggregated ROC comparison plot for ' label '...']);
        plot_aggregated_roc_comparison(current_data, label);

        disp(['Generating aggregated SVM classification plot for ' label '...']);
        plot_aggregated_svm(current_data, label);

        disp(['Generating aggregated ANOVA plot for ' label '...']);
        plot_aggregated_anova(current_data, label);

        disp(['Generating aggregated baseline comparison plot for ' label '...']);
        plot_aggregated_baseline_comparison(current_data, label);
    else
        warning('Data for type ''%s'' not found in aggregated file. Skipping.', type);
    end
end

% End timer and provide user feedback
toc;
disp('Plotting pipeline finished.');
