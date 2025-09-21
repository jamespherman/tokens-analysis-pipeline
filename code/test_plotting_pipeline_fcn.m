function test_plotting_pipeline_fcn(unique_id)
%% test_plotting_pipeline_fcn.m
%
%   This function tests the aggregated plotting pipeline. It verifies that
%   the main plotting scripts run without errors.
%
%   Args:
%       unique_id (char): The unique identifier for the session to be tested.
%                         Note: This is used to load the correct manifest, but
%                         the plotting scripts themselves run on aggregated
%                         data.
%
% Author: Jules
% Date: 2025-09-21

% Start timer
tic;

% In-line function to report timing
giveFeed = @(x)disp([num2str(toc) 's - ' x]);

%% Setup Paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Intro
giveFeed(sprintf('Testing plotting pipeline for session: %s', unique_id));

%% Load Manifest
giveFeed('Step 1: Loading session manifest...');
project_root = fullfile(findOneDrive, 'Code', 'tokens-analysis-pipeline');
manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');
manifest = readtable(manifest_path);

% Find the row for the current session
session_idx = find(strcmp(manifest.unique_id, unique_id));
if isempty(session_idx)
    error('test_plotting_pipeline:sessionNotFound', ...
          'Session with unique_id "%s" not found in the manifest.', ...
          unique_id);
end

%% Run Plotting Pipeline
% This function will execute the aggregated plotting scripts.
% It assumes that the analysis pipeline has been run and aggregated data
% files are available in `data/processed`.

% Step 1: Run run_plotting_pipeline.m
try
    giveFeed('Running run_plotting_pipeline...');
    run_plotting_pipeline(unique_id);
    giveFeed('run_plotting_pipeline completed.');
catch ME
    fprintf(2, 'Error running run_plotting_pipeline for %s\n', unique_id);
    rethrow(ME);
end

%% Done
giveFeed('Plotting pipeline test complete.');

end
