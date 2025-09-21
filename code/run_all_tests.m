%% run_all_tests.m
%
%   This script serves as the master test runner for the entire tokens-analysis-pipeline.
%   It iterates through a predefined list of sessions and calls a series of
%   test functions to verify the functionality of each major component of the
%   pipeline, from data preparation to final plotting.
%
%   Usage:
%       - Ensure that the MATLAB path is set to include the 'code' and 'code/utils'
%         directories.
%       - Run this script from the MATLAB command window.
%
%   Tests Executed:
%       1. test_neuron_diagnostics_fcn: Verifies neuron screening and diagnostic PDF generation.
%       2. test_core_data_preparation_fcn: Verifies the core data preparation pipeline.
%       3. test_analysis_pipeline_fcn: Verifies the main analysis functions.
%       4. test_plotting_pipeline_fcn: Verifies the aggregated plotting functions.
%
% Author: Jules
% Date: 2025-09-21

%% Setup
clear; clc; close all;

% Add utility functions to the MATLAB path
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir));
addpath(fullfile(script_dir, 'utils'));

%% Define Test Sessions
% List of unique_ids for the sessions to be tested.
% These sessions are chosen to cover different conditions (e.g., SNc, SC).
test_sessions = { ...
    'Feynman_08_15_2025_SC', ...
    'Feynman_08_11_2025_SNc' ...
};

%% Run Tests
fprintf('Starting test suite...\n');
fprintf('--------------------------------\n');

for i_session = 1:numel(test_sessions)
    unique_id = test_sessions{i_session};
    fprintf('Running tests for session: %s\n', unique_id);

    try
        % Test 1: Neuron Diagnostics
        fprintf('  [1/4] Running neuron diagnostics test...\n');
        test_neuron_diagnostics_fcn(unique_id);

        % Test 2: Core Data Preparation
        fprintf('  [2/4] Running core data preparation test...\n');
        test_core_data_preparation_fcn(unique_id);

        % Test 3: Analysis Pipeline
        fprintf('  [3/4] Running analysis pipeline test...\n');
        test_analysis_pipeline_fcn(unique_id);

        % Test 4: Plotting Pipeline
        fprintf('  [4/4] Running plotting pipeline test...\n');
        test_plotting_pipeline_fcn(unique_id);

        fprintf('Tests for %s completed successfully.\n', unique_id);
        fprintf('--------------------------------\n');

    catch ME
        fprintf(2, 'Error running tests for session: %s\n', unique_id);
        fprintf(2, 'Error message: %s\n', ME.message);
        fprintf(2, 'Stack trace:\n');
        for i_stack = 1:numel(ME.stack)
            fprintf(2, '  File: %s, Name: %s, Line: %d\n', ...
                ME.stack(i_stack).file, ME.stack(i_stack).name, ME.stack(i_stack).line);
        end
        fprintf('--------------------------------\n');
    end
end

fprintf('Test suite finished.\n');
