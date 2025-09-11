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

% --- USER TOGGLES ---
% Force re-computation of all steps, even if marked 'complete'
force_rerun = false;
% --- END USER TOGGLES ---

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

%% Load Analysis Plan
giveFeed('Loading analysis plan...');
analysis_plan = define_analysis_plan();
giveFeed('Analysis plan loaded.');

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

    % >> SESSION-LEVEL PROGRESS REPORTING (PART B)
    % --- Dry run to calculate total number of steps for this session ---
    n_total_steps = 0;
    unique_id = session_id; % Use a consistent variable name for fprintf

    % 1. Screening
    if ~strcmp(manifest.screening_status{i}, 'complete') || force_rerun
        n_total_steps = n_total_steps + 1;
    end

    % 2. PDF Generation
    diag_output_dir_dry_run = fullfile(project_root, 'figures', unique_id);
    if (~exist(diag_output_dir_dry_run, 'dir') || isempty(dir(fullfile(diag_output_dir_dry_run, '*.pdf')))) || force_rerun
        n_total_steps = n_total_steps + 1;
    end

    % 3. Data Prep
    if ~strcmp(manifest.dataprep_status{i}, 'complete') || force_rerun
        n_total_steps = n_total_steps + 1;
    end

    % 4. Analyses
    analysis_name_baseline = 'baseline_comparison';
    if isfield(analysis_plan, analysis_name_baseline)
        conditions_to_run = analysis_plan.(analysis_name_baseline).conditions_to_run;
        for j = 1:length(conditions_to_run)
            condition_name = conditions_to_run{j};
            if (~isfield(session_data, 'analysis') || ...
               ~isfield(session_data.analysis, analysis_name_baseline) || ...
               ~isfield(session_data.analysis.(analysis_name_baseline), condition_name)) || force_rerun
                n_total_steps = n_total_steps + 1;
            end
        end
    end
    analysis_name_roc = 'roc_comparison';
    if isfield(analysis_plan, analysis_name_roc)
        comparisons_to_run = analysis_plan.(analysis_name_roc).comparisons_to_run;
        for j = 1:length(comparisons_to_run)
            comp = comparisons_to_run(j);
            if (~isfield(session_data, 'analysis') || ...
               ~isfield(session_data.analysis, analysis_name_roc) || ...
               ~isfield(session_data.analysis.(analysis_name_roc), comp.name)) || force_rerun
                n_total_steps = n_total_steps + 1;
            end
        end
    end

    step_counter = 0;
    % --- End of dry run ---

    % --- 1. Neuron Screening ---
    if ~strcmp(manifest.screening_status{i}, 'complete') || force_rerun
        step_counter = step_counter + 1;
        fprintf('\n--- Session %s: Starting Step %d of %d: Neuron Screening ---\n', unique_id, step_counter, n_total_steps);
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

    % --- Per-Neuron Diagnostic PDF Generation ---
    giveFeed('Checking for per-neuron diagnostic PDFs...');
    diag_output_dir = fullfile(project_root, 'figures', session_id);

    % Check if the directory exists and contains any PDF files
    if (~exist(diag_output_dir, 'dir') || isempty(dir(fullfile(diag_output_dir, '*.pdf')))) || force_rerun
        step_counter = step_counter + 1;
        fprintf('\n--- Session %s: Starting Step %d of %d: Diagnostic PDF Generation ---\n', unique_id, step_counter, n_total_steps);
        giveFeed('Generating diagnostic PDF...');
        if ~exist(diag_output_dir, 'dir')
            mkdir(diag_output_dir);
        end
        generate_neuron_summary_pdf(session_data, selected_neurons, ...
            session_id, diag_output_dir);
        giveFeed('Diagnostic PDF generation complete.');
    else
        giveFeed('Diagnostic PDF already exists, skipping...');
    end

    % --- 2. Core Data Preparation ---
    if ~strcmp(manifest.dataprep_status{i}, 'complete') || force_rerun
        step_counter = step_counter + 1;
        fprintf('\n--- Session %s: Starting Step %d of %d: Core Data Preparation ---\n', unique_id, step_counter, n_total_steps);
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

    % --- 4. On-Demand Analysis Execution ---
    giveFeed('Checking for missing analyses...');
    data_updated = false; % Flag to track if we modify session_data
    all_analyses_complete = true; % Flag to track if all analyses are present

    % A. Baseline Comparison Analyses
    analysis_name = 'baseline_comparison';
    if isfield(analysis_plan, analysis_name)
        conditions_to_run = analysis_plan.(analysis_name).conditions_to_run;
        for j = 1:length(conditions_to_run)
            condition_name = conditions_to_run{j};
            if (~isfield(session_data, 'analysis') || ...
               ~isfield(session_data.analysis, analysis_name) || ...
               ~isfield(session_data.analysis.(analysis_name), condition_name)) || force_rerun

                step_counter = step_counter + 1;
                progress_message = sprintf('Baseline Comparison for %s', condition_name);
                fprintf('\n--- Session %s: Starting Step %d of %d: %s ---\n', unique_id, step_counter, n_total_steps, progress_message);

                all_analyses_complete = false; % Mark as incomplete
                giveFeed(sprintf('--> Running missing analysis: %s for %s', ...
                    analysis_name, condition_name));

                analysis_result = analyze_baseline_comparison(core_data, ...
                    conditions, is_av_session, 'condition', condition_name);

                session_data.analysis.(analysis_name).(condition_name) = analysis_result;
                data_updated = true;
            end
        end
    end

    % B. ROC Comparison Analyses
    analysis_name = 'roc_comparison';
    if isfield(analysis_plan, analysis_name)
        comparisons_to_run = analysis_plan.(analysis_name).comparisons_to_run;
        for j = 1:length(comparisons_to_run)
            comp = comparisons_to_run(j);
            if (~isfield(session_data, 'analysis') || ...
               ~isfield(session_data.analysis, analysis_name) || ...
               ~isfield(session_data.analysis.(analysis_name), comp.name)) || force_rerun

                step_counter = step_counter + 1;
                progress_message = sprintf('ROC Comparison for %s', comp.name);
                fprintf('\n--- Session %s: Starting Step %d of %d: %s ---\n', unique_id, step_counter, n_total_steps, progress_message);

                all_analyses_complete = false; % Mark as incomplete
                giveFeed(sprintf('--> Running missing analysis: %s for %s', ...
                    analysis_name, comp.name));

                analysis_result = analyze_roc_comparison(core_data, ...
                    conditions, is_av_session, 'comparison', comp);

                session_data.analysis.(analysis_name).(comp.name) = analysis_result;
                data_updated = true;
            end
        end
    end

    % --- 5. Save Updated Data ---
    if data_updated
        giveFeed('Analyses run, saving updated session_data.mat...');
        save(session_data_path, 'session_data', '-v7.3');
        giveFeed('Save complete.');
    else
        giveFeed('No new analyses were required for this session.');
    end

    % Update manifest status only if all analyses are confirmed complete
    if all_analyses_complete
        manifest.analysis_status{i} = 'complete';
    end

    giveFeed(sprintf('--- Finished processing for session: %s ---\n', session_id));
end

%% Finalize
giveFeed('All sessions processed. Saving updated manifest...');
writetable(manifest, manifest_path);
giveFeed('Manifest saved. Total time elapsed:');
toc
