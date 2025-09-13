%% run_tokens_analysis.m
%
% Main script to run the entire tokens task analysis pipeline. It iterates
% through a session manifest, performs neuron screening, prepares core
% data, and runs a series of specified analyses based on a centralized plan.
%
% The script is designed to be idempotent; it checks the status of each
% step in the manifest and skips steps that are already marked 'complete'.
%
% Author: Jules
% Date: 2025-09-12 (Refactored)
%

%% Setup
clear; clc; close all;

% --- USER TOGGLES ---
force_rerun = struct(...
    'screening', false, ...
    'diag_pdfs', false, ...
    'dataprep',  false, ...
    'analyses',  false ...
);
% --- END USER TOGGLES ---

[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);
addpath(project_root);

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
% The single source of truth for the analysis plan is now
% define_task_conditions. We call it without arguments to get the plan.
[~, ~, analysis_plan] = define_task_conditions();
giveFeed('Analysis plan loaded.');

%% Iterate Through Sessions
for i = 1:height(manifest)
    session_id = manifest.unique_id{i};
    giveFeed(sprintf('--- Starting processing for session: %s ---', ...
        session_id));

    one_drive_path = findOneDrive;
    session_data_path = fullfile(one_drive_path, ...
        'Neuronal Data Analysis', session_id, [session_id ...
        '_session_data.mat']);

    if ~exist(session_data_path, 'file')
        warning('run_tokens_analysis:sessionDataNotFound', ...
            'session_data.mat not found for %s. Skipping.', ...
            session_id);
        continue;
    end

    giveFeed(sprintf('Loading data for %s...', session_id));
    load(session_data_path, 'session_data');
    giveFeed('Data loaded.');

    data_updated = false;

    % --- Define Task Conditions for this session ---
    giveFeed('Defining task conditions...');
    [conditions, is_av_session] = define_task_conditions( ...
        session_data.trialInfo, ...
        session_data.eventTimes, session_data.metadata.unique_id);
    giveFeed('Task conditions defined.');

    % --- Dry run to calculate total number of steps for this session ---
    n_total_steps = 0;
    unique_id = session_id;
    if ~isfield(session_data, 'metrics') || force_rerun.screening
        n_total_steps = n_total_steps + 1;
    end
    if ~strcmp(manifest.screening_status{i}, 'complete') || ...
            force_rerun.screening
        n_total_steps = n_total_steps + 1;
    end
    diag_output_dir_dry_run = fullfile(project_root, 'figures', unique_id);
    if (~exist(diag_output_dir_dry_run, 'dir') || ...
            isempty(dir(fullfile(diag_output_dir_dry_run, '*.pdf'))))
        n_total_steps = n_total_steps + 1;
    end
    if ~strcmp(manifest.dataprep_status{i}, 'complete') || ...
            force_rerun.dataprep
        n_total_steps = n_total_steps + 1;
    end

    % Count planned analyses
    for j = 1:length(analysis_plan.baseline_plan)
        comp = analysis_plan.baseline_plan(j);
        if (~isfield(session_data, 'analysis') || ...
           ~isfield(session_data.analysis, 'baseline_comparison') || ...
           ~isfield(session_data.analysis.baseline_comparison, ...
           comp.name)) || force_rerun.analyses
            n_total_steps = n_total_steps + 1;
        end
    end
    for j = 1:length(analysis_plan.roc_plan)
        comp = analysis_plan.roc_plan(j);
        if (~isfield(session_data, 'analysis') || ...
           ~isfield(session_data.analysis, 'roc_comparison') || ...
           ~isfield(session_data.analysis.roc_comparison, comp.name)) ...
           || force_rerun.analyses
            n_total_steps = n_total_steps + 1;
        end
    end
    if analysis_plan.anova_plan.run
        if (~isfield(session_data, 'analysis') || ...
           ~isfield(session_data.analysis, 'anova_results')) || ...
            force_rerun.analyses
            n_total_steps = n_total_steps + 1;
        end
    end
    step_counter = 0;

    % --- Pipeline Stages ---
    % ... (screening, PDF gen, data prep stages are unchanged) ...
        % --- 1. Neuron Screening ---
    if ~strcmp(manifest.screening_status{i}, 'complete') || ...
        force_rerun.screening
        step_counter = step_counter + 1;
        fprintf(['\n--- Session %s: Starting Step %d of %d: ' ...
            '        Neuron Screening ---\n'], unique_id, ...
            step_counter, n_total_steps);
        giveFeed('Screening status is ''pending''. Running screening...');

        session_data.metadata.unique_id = session_id;
        if contains(session_id, 'SNc')
            selected_neurons = screen_da_neurons(session_data, session_id);
        elseif contains(session_id, 'SC')
            [selected_neurons, sig_epoch_comp, scSide] = ...
                screen_sc_neurons(session_data);
            session_data.analysis.scSide = scSide;
            session_data.analysis.sig_epoch_comparison = sig_epoch_comp;
        else
            warning('run_tokens_analysis:unknownSessionType', ...
                'Unknown session type for %s. Cannot screen neurons.', ...
                session_id);
            continue;
        end

        session_data.analysis.selected_neurons = selected_neurons;

        data_updated = true; % Mark data as updated
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
    if (~exist(diag_output_dir, 'dir') || isempty(dir(fullfile( ...
            diag_output_dir, '*.pdf')))) || force_rerun.diag_pdfs 
        step_counter = step_counter + 1;
        fprintf(['\n--- Session %s: Starting Step %d of %d: ' ...
            'Diagnostic PDF Generation ---\n'], unique_id, ...
            step_counter, n_total_steps);
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
    if ~strcmp(manifest.dataprep_status{i}, 'complete') || ...
        force_rerun.dataprep
        step_counter = step_counter + 1;
        fprintf(['\n--- Session %s: Starting Step %d of %d: ' ...
            'Core Data Preparation ---\n'], unique_id, ...
            step_counter, n_total_steps);
        giveFeed(['Data prep status is ''pending''. ' ...
            '            Running prepare_core_data...']);
        core_data = prepare_core_data(session_data, selected_neurons);
        session_data.analysis.core_data = core_data;

        data_updated = true; % Mark data as updated
        manifest.dataprep_status{i} = 'complete';
        giveFeed('Data prep complete.');
    else
        giveFeed('Data prep already complete. Loading core_data.');
        core_data = session_data.analysis.core_data;
    end

    % --- On-Demand Analysis Execution ---
    giveFeed('Checking for missing analyses...');

    % A. Baseline Comparison Analyses
    for j = 1:length(analysis_plan.baseline_plan)
        comp = analysis_plan.baseline_plan(j);
        if comp.is_av_only && ~is_av_session, continue; end

        if (~isfield(session_data, 'analysis') || ...
           ~isfield(session_data.analysis, 'baseline_comparison') || ...
           ~isfield(session_data.analysis.baseline_comparison, ...
           comp.name)) || force_rerun.analyses
            step_counter = step_counter + 1;
            fprintf(['\n--- Session %s: Step %d/%d: Baseline ' ...
                'Comparison for %s ---\n'], unique_id, step_counter, ...
                n_total_steps, comp.name);
            giveFeed(sprintf('--> Running Baseline Comparison: %s', ...
                comp.name));
            result = analyze_baseline_comparison(core_data, conditions, ...
                is_av_session, 'condition', comp.name);
            session_data.analysis.baseline_comparison.(comp.name) = result;
            data_updated = true;
        end
    end

    % B. ROC Comparison Analyses
    for j = 1:length(analysis_plan.roc_plan)
        comp = analysis_plan.roc_plan(j);
        if comp.is_av_only && ~is_av_session, continue; end

        if (~isfield(session_data, 'analysis') || ...
           ~isfield(session_data.analysis, 'roc_comparison') || ...
           ~isfield(session_data.analysis.roc_comparison, comp.name)) ...
           || force_rerun.analyses
            step_counter = step_counter + 1;
            fprintf(['\n--- Session %s: Step %d/%d: ROC Comparison ' ...
                'for %s ---\n'], unique_id, step_counter, ...
                n_total_steps, comp.name);
            giveFeed(sprintf('--> Running ROC Comparison: %s', comp.name));
            result = analyze_roc_comparison(core_data, conditions, ...
                is_av_session, 'comparison', comp);
            session_data.analysis.roc_comparison.(comp.name) = ...
                result.(comp.name);
            data_updated = true;
        end
    end

    % C. N-way ANOVA Analysis
    if analysis_plan.anova_plan.run
        if (~isfield(session_data, 'analysis') || ...
           ~isfield(session_data.analysis, 'anova_results')) ...
           || force_rerun.analyses
            step_counter = step_counter + 1;
            fprintf('\n--- Session %s: Step %d/%d: N-way ANOVA ---\n', ...
                unique_id, step_counter, n_total_steps);
            giveFeed('--> Running N-way ANOVA');
            session_data = analyze_anova(session_data, core_data, ...
                conditions);
            data_updated = true;
        end
    end

    % --- Save Updated Data ---
    if data_updated
        giveFeed('Data was updated, saving back to session_data.mat...');
        save(session_data_path, 'session_data', '-v7.3');
        giveFeed('Save complete.');
    else
        giveFeed(['No new analyses or processing were required for ' ...
            '                this session.']);
    end

    % --- Verify Analysis Completion & Update Manifest ---
    is_analysis_complete = true;
    for j = 1:length(analysis_plan.baseline_plan)
        comp = analysis_plan.baseline_plan(j);
        if ~isfield(session_data, 'analysis') || ...
           ~isfield(session_data.analysis, 'baseline_comparison') || ...
           ~isfield(session_data.analysis.baseline_comparison, comp.name)
            is_analysis_complete = false; break;
        end
    end
    if is_analysis_complete
        for j = 1:length(analysis_plan.roc_plan)
            comp = analysis_plan.roc_plan(j);
            if ~isfield(session_data, 'analysis') || ...
               ~isfield(session_data.analysis, 'roc_comparison') || ...
               ~isfield(session_data.analysis.roc_comparison, comp.name)
                is_analysis_complete = false; break;
            end
        end
    end
    if is_analysis_complete && analysis_plan.anova_plan.run
        if ~isfield(session_data, 'analysis') || ...
           ~isfield(session_data.analysis, 'anova_results')
            is_analysis_complete = false;
        end
    end

    if is_analysis_complete
        manifest.analysis_status{i} = 'complete';
    end

    giveFeed(sprintf('--- Finished processing for session: %s ---\n', ...
        session_id));
end

%% Finalize
giveFeed('All sessions processed. Saving updated manifest...');
writetable(manifest, manifest_path);
giveFeed('Manifest saved. Total time elapsed:');
toc
