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
load(aggFileName);
disp('Data loaded successfully.');

% Call plotting functions for SC data
disp('--- Generating SC Plots ---');
disp('Generating aggregated ROC comparison plot for SC...');
plot_aggregated_roc_comparison(aggregated_sc_data, 'SC');
disp('Generating aggregated ANOVA plot for SC...');
plot_aggregated_anova(aggregated_sc_data, 'SC');
disp('Generating aggregated baseline comparison plot for SC...');
plot_aggregated_baseline_comparison(aggregated_sc_data, 'SC');

% Call plotting functions for SNc data
disp('--- Generating SNc Plots ---');
disp('Generating aggregated ROC comparison plot for SNc...');
plot_aggregated_roc_comparison(aggregated_snc_data, 'SNc');
disp('Generating aggregated ANOVA plot for SNc...');
plot_aggregated_anova(aggregated_snc_data, 'SNc');
disp('Generating aggregated baseline comparison plot for SNc...');
plot_aggregated_baseline_comparison(aggregated_snc_data, 'SNc');

% End timer and provide user feedback
toc;
disp('Plotting pipeline finished.');
