function test_analysis_pipeline_fcn(unique_id)
%% test_analysis_pipeline_fcn.m
%
%   This function tests the core analysis pipeline for a single session. It
%   verifies that the main analysis scripts run without errors.
%
%   Args:
%       unique_id (char): The unique identifier for the session to be tested.
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
giveFeed(sprintf('Testing analysis pipeline for session: %s', unique_id));

%% Load Manifest
giveFeed('Step 1: Loading session manifest...');
project_root = fullfile(findOneDrive, 'Code', 'tokens-analysis-pipeline');
manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');
manifest = readtable(manifest_path);

% Find the row for the current session
session_idx = find(strcmp(manifest.unique_id, unique_id));
if isempty(session_idx)
    error('test_analysis_pipeline:sessionNotFound', ...
          'Session with unique_id "%s" not found in the manifest.', ...
          unique_id);
end

%% Run Analysis Pipeline
% This function will execute the core analysis scripts in sequence.
% It assumes that the raw data has been processed and `session_data.mat` and
% `core_data` are available.

% Step 1: Run run_tokens_analysis.m
try
    giveFeed('Running run_tokens_analysis...');
    run_tokens_analysis(unique_id);
    giveFeed('run_tokens_analysis completed.');
catch ME
    fprintf(2, 'Error running run_tokens_analysis for %s\n', unique_id);
    rethrow(ME);
end

% Step 2: Run analyze_baseline_comparison.m
try
    giveFeed('Running analyze_baseline_comparison...');
    analyze_baseline_comparison(unique_id);
    giveFeed('analyze_baseline_comparison completed.');
catch ME
    fprintf(2, 'Error running analyze_baseline_comparison for %s\n', unique_id);
    rethrow(ME);
end

% Step 3: Run analyze_roc_comparison.m
try
    giveFeed('Running analyze_roc_comparison...');
    analyze_roc_comparison(unique_id);
    giveFeed('analyze_roc_comparison completed.');
catch ME
    fprintf(2, 'Error running analyze_roc_comparison for %s\n', unique_id);
    rethrow(ME);
end

% Step 4: Run analyze_anova.m
try
    giveFeed('Running analyze_anova...');
    analyze_anova(unique_id);
    giveFeed('analyze_anova completed.');
catch ME
    fprintf(2, 'Error running analyze_anova for %s\n', unique_id);
    rethrow(ME);
end

%% Done
giveFeed('Analysis pipeline test complete.');

end
