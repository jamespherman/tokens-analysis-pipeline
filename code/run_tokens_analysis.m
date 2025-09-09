%% run_tokens_analysis.m
%
% Main script to run the entire tokens task analysis pipeline. It iterates
% through a session manifest, performs neuron screening, prepares core
% data, and runs a series of specified analyses.
%
% The script is designed to be idempotent; it checks the status of each
% step in the manifest and skips steps that are already marked 'complete'.
%
% Author: Jules
% Date: 2025-09-09
%

%% Setup
clear; clc; close all;

% Add utility functions to the MATLAB path
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir); % Assumes script is in code/
addpath(project_root); % Add project root to path

% Start timer and provide feedback
tic;
giveFeed = @(x)disp([num2str(round(toc, 1)) 's - ' x]);

%% Load Manifest
giveFeed('Loading session manifest...');
manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');
if ~exist(manifest_path, 'file')
    error('run_tokens_analysis:manifestNotFound', ...
        'Session manifest not found at: %s', manifest_path);
end
manifest = readtable(manifest_path);
giveFeed('Manifest loaded.');

%% Iterate Through Sessions
for i = 1:height(manifest)
    session_id = manifest.unique_id{i};
    giveFeed(sprintf('--- Starting processing for session: %s ---', session_id));

    % Define path to the session_data.mat file
    one_drive_path = findOneDrive;
    session_data_path = fullfile(one_drive_path, ...
        'Neuronal Data Analysis', session_id, ...
        [session_id '_session_data.mat']);

    if ~exist(session_data_path, 'file')
        warning('run_tokens_analysis:sessionDataNotFound', ...
            'session_data.mat not found for %s. Skipping.', session_id);
        continue;
    end

    % Load the session data
    giveFeed(sprintf('Loading data for %s...', session_id));
    load(session_data_path, 'session_data');
    giveFeed('Data loaded.');

    % --- 1. Neuron Screening ---
    if ~strcmp(manifest.screening_status{i}, 'complete')
        giveFeed('Screening status is ''pending''. Running screening...');

        session_data.metadata.unique_id = session_id;
        if contains(session_id, 'SNc')
            selected_neurons = screen_da_neurons(session_data, session_id);
        elseif contains(session_id, 'SC')
            [selected_neurons, sig_epoch_comp, scSide] = screen_sc_neurons(session_data);
            session_data.analysis.scSide = scSide;
            session_data.analysis.sig_epoch_comparison = sig_epoch_comp;
        else
            warning('run_tokens_analysis:unknownSessionType', ...
                'Unknown session type for %s. Cannot screen neurons.', session_id);
            continue;
        end

        session_data.analysis.selected_neurons = selected_neurons;

        giveFeed('Saving screening results back to session_data.mat...');
        save(session_data_path, 'session_data', '-v7.3');
        manifest.screening_status{i} = 'complete';
        giveFeed('Screening complete.');
    else
        giveFeed('Screening already complete. Loading results.');
        selected_neurons = session_data.analysis.selected_neurons;
    end

    % --- 2. Core Data Preparation ---
    if ~strcmp(manifest.dataprep_status{i}, 'complete')
        giveFeed('Data prep status is ''pending''. Running prepare_core_data...');
        core_data = prepare_core_data(session_data, selected_neurons);
        session_data.analysis.core_data = core_data;

        giveFeed('Saving core_data back to session_data.mat...');
        save(session_data_path, 'session_data', '-v7.3');
        manifest.dataprep_status{i} = 'complete';
        giveFeed('Data prep complete.');
    else
        giveFeed('Data prep already complete. Loading core_data.');
        core_data = session_data.analysis.core_data;
    end

    % --- 3. Define Task Conditions ---
    giveFeed('Defining task conditions...');
    [conditions, is_av_session] = define_task_conditions(session_data.trialInfo, ...
        session_data.eventTimes, session_data.metadata.unique_id);
    giveFeed('Task conditions defined.');

    % --- 4. Run Analyses ---
    if ~isfield(manifest, 'analysis_status')
        manifest.analysis_status = repmat({'pending'}, height(manifest), 1);
    end

    if ~strcmp(manifest.analysis_status{i}, 'complete')
        giveFeed('Analysis status is ''pending''. Running analyses...');

        % Run Baseline vs. Post-Event Analysis
        giveFeed('--> Running baseline comparison analysis...');
        results_baseline = analyze_baseline_comparison(core_data, conditions, is_av_session);

        % Run Between-Condition ROC Analysis
        giveFeed('--> Running ROC comparison analysis...');
        results_roc = analyze_roc_comparison(core_data, conditions, is_av_session);

        % --- 5. Save Analysis Results ---
        output_dir = fullfile(project_root, 'data', 'processed', session_id);
        if ~exist(output_dir, 'dir')
           mkdir(output_dir);
           giveFeed(sprintf('Created output directory: %s', output_dir));
        end

        results_path = fullfile(output_dir, 'analysis_results.mat');
        giveFeed(sprintf('Saving analysis results to: %s', results_path));
        save(results_path, 'results_baseline', 'results_roc', 'conditions');

        manifest.analysis_status{i} = 'complete';
        giveFeed('Analyses complete and results saved.');
    else
        giveFeed('Analyses already complete. Skipping.');
    end

    % --- 6. Update Manifest File ---
    giveFeed(sprintf('Updating manifest for %s...', session_id));
    writetable(manifest, manifest_path);
    giveFeed(sprintf('--- Finished processing for session: %s ---\n', session_id));
end

giveFeed('All sessions processed. Total time elapsed:');
toc
