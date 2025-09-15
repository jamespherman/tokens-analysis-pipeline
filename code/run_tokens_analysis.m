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
    'screening', true, ...
    'diag_pdfs', false, ...
    'dataprep',  true, ...
    'analyses',  true ...
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

%% Dynamically Generate Required Alignment Events
giveFeed('Generating required alignment events list from analysis plan...');
all_events = {};
plan_fields = fieldnames(analysis_plan);
for i = 1:length(plan_fields)
    sub_plan = analysis_plan.(plan_fields{i});
    if isstruct(sub_plan)
        for j = 1:length(sub_plan)
            if isfield(sub_plan(j), 'event')
                all_events{end+1} = sub_plan(j).event;
            end
        end
    end
end
alignment_events = unique(all_events);
giveFeed('Alignment events list generated.');

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
        path_to_check = fullfile('analysis', 'baseline_comparison', 'CUE_ON', comp.name);
        S = substruct('.', strsplit(path_to_check, '/'));
        try
            subsref(session_data, S);
        catch
            n_total_steps = n_total_steps + 1;
        end
    end
    for j = 1:length(analysis_plan.roc_plan)
        comp = analysis_plan.roc_plan(j);
        path_to_check = fullfile('analysis', 'roc_comparison', comp.event, comp.name);
        S = substruct('.', strsplit(path_to_check, '/'));
        try
            subsref(session_data, S);
        catch
            n_total_steps = n_total_steps + 1;
        end
    end
    anova_step_counted = false;
    for j = 1:length(analysis_plan.anova_plan)
        plan = analysis_plan.anova_plan(j);
        if plan.run && ~anova_step_counted
            if (~isfield(session_data, 'analysis') || ...
               ~isfield(session_data.analysis, 'anova_results')) || ...
                force_rerun.analyses
                n_total_steps = n_total_steps + 1;
                anova_step_counted = true; % Count only once
            end
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
        core_data = prepare_core_data(session_data, selected_neurons, alignment_events);
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
    % Iterate through the plan, running analysis for each condition.
    for j = 1:length(analysis_plan.baseline_plan)
        comp_plan = analysis_plan.baseline_plan(j);
        if comp_plan.is_av_only && ~is_av_session, continue; end

        % Check if this specific analysis needs to be run.
        % We check for the result of the first alignment event ('CUE_ON')
        % as a proxy for the entire analysis set for this condition.
        analysis_needed = force_rerun.analyses;
        if ~analysis_needed
            try
                % Check for existence of CUE_ON result for this condition
                subsref(session_data, substruct('.', 'analysis', ...
                    '.', 'baseline_comparison', '.', 'CUE_ON', ...
                    '.', comp_plan.name));
            catch
                analysis_needed = true;
            end
        end

        if analysis_needed
            step_counter = step_counter + 1;
            fprintf(['\n--- Session %s: Step %d/%d: Baseline ' ...
                'Comparison for %s ---\n'], unique_id, step_counter, ...
                n_total_steps, comp_plan.name);
            giveFeed(sprintf(['--> Running analyze_baseline_comparison ' ...
                'for: %s'], comp_plan.name));

            % This function returns a struct with event names as fields.
            result_by_event = analyze_baseline_comparison(core_data, ...
                conditions, is_av_session, 'condition', comp_plan.name);

            % Merge the results into the standardized structure.
            event_names = fieldnames(result_by_event);
            for k = 1:length(event_names)
                event_name = event_names{k};
                % Correct the assignment to prevent redundant nesting.
                % The result from the analysis function is already nested
                % by event and condition, so we need to extract the
                % specific condition's data.
                session_data.analysis.baseline_comparison.(event_name).(comp_plan.name) = ...
                    result_by_event.(event_name).(comp_plan.name);
            end
            data_updated = true;
        end
    end

    % B. ROC Comparison Analyses
    % Iterate through each comparison defined in the plan.
    for j = 1:length(analysis_plan.roc_plan)
        comp_plan = analysis_plan.roc_plan(j);
        if comp_plan.is_av_only && ~is_av_session, continue; end

        % Check if this specific analysis needs to be run.
        analysis_needed = force_rerun.analyses;
        if ~analysis_needed
            try
                subsref(session_data, substruct('.', 'analysis', ...
                    '.', 'roc_comparison', '.', comp_plan.event, ...
                    '.', comp_plan.name));
            catch
                analysis_needed = true;
            end
        end

        if analysis_needed
            step_counter = step_counter + 1;
            fprintf(['\n--- Session %s: Step %d/%d: ROC Comparison ' ...
                'for %s ---\n'], unique_id, step_counter, ...
                n_total_steps, comp_plan.name);
            giveFeed(sprintf('--> Running analyze_roc_comparison for: %s', ...
                comp_plan.name));

            % The function returns a result nested by event and comp name.
            result = analyze_roc_comparison(core_data, conditions, ...
                is_av_session, 'comparison', comp_plan);

            % Extract the specific result and store it.
            session_data.analysis.roc_comparison.(comp_plan.event).(comp_plan.name) = ...
                result.(comp_plan.event).(comp_plan.name);
            data_updated = true;
        end
    end

    % C. N-way ANOVA Analysis
    % This analysis runs once for all events defined internally.
    % Check if it's needed.
    analysis_needed = force_rerun.analyses;
    if ~analysis_needed
        try
            % Check for the existence of the top-level anova_results field.
            subsref(session_data, substruct('.', 'analysis', ...
                '.', 'anova_results'));
        catch
            analysis_needed = true;
        end
    end

    if analysis_needed
        step_counter = step_counter + 1;
        fprintf('\n--- Session %s: Step %d/%d: N-way ANOVA ---\n', ...
            unique_id, step_counter, n_total_steps);
        giveFeed('--> Running analyze_anova...');

        % The function takes session_data and returns it with results added.
        session_data = analyze_anova(session_data, core_data, conditions);
        data_updated = true;
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
    % This block verifies that all analyses defined in the analysis_plan
    % have been successfully run and their results saved to the
    % session_data struct.
    is_analysis_complete = true; % Assume complete until proven otherwise

    % 1. Verify Baseline Comparison Results
    for j = 1:length(analysis_plan.baseline_plan)
        plan = analysis_plan.baseline_plan(j);
        if plan.is_av_only && ~is_av_session, continue; end

        % This analysis runs for all alignment events. We must check all.
        for k = 1:length(alignment_events)
            event_name = alignment_events{k};
            try
                subsref(session_data, substruct('.', 'analysis', ...
                    '.', 'baseline_comparison', '.', event_name, ...
                    '.', plan.name));
            catch
                is_analysis_complete = false;
                fprintf('Manifest Check Failed: Missing Baseline Comparison -> Event: %s, Cond: %s\n', ...
                    event_name, plan.name);
                break; % Exit inner loop
            end
        end
        if ~is_analysis_complete, break; end % Exit outer loop
    end

    % 2. Verify ROC Comparison Results
    if is_analysis_complete
        for j = 1:length(analysis_plan.roc_plan)
            plan = analysis_plan.roc_plan(j);
            if plan.is_av_only && ~is_av_session, continue; end

            try
                subsref(session_data, substruct('.', 'analysis', ...
                    '.', 'roc_comparison', '.', plan.event, ...
                    '.', plan.name));
            catch
                is_analysis_complete = false;
                fprintf('Manifest Check Failed: Missing ROC Comparison -> Event: %s, Comp: %s\n', ...
                    plan.event, plan.name);
                break;
            end
        end
    end

    % 3. Verify N-way ANOVA Results
    if is_analysis_complete
        for j = 1:length(analysis_plan.anova_plan)
            plan = analysis_plan.anova_plan(j);
            if ~plan.run, continue; end

            try
                subsref(session_data, substruct('.', 'analysis', ...
                    '.', 'anova_results', '.', plan.event));
            catch
                is_analysis_complete = false;
                fprintf('Manifest Check Failed: Missing ANOVA -> Event: %s\n', ...
                    plan.event);
                break;
            end
        end
    end

    % Update manifest status based on verification outcome
    if is_analysis_complete
        manifest.analysis_status{i} = 'complete';
        giveFeed('All planned analyses are complete.');
    else
        manifest.analysis_status{i} = 'pending';
        giveFeed('Verification failed. Some analyses are missing.');
    end

    giveFeed(sprintf('--- Finished processing for session: %s ---\n', ...
        session_id));
end

%% Finalize
giveFeed('All sessions processed. Saving updated manifest...');
writetable(manifest, manifest_path);
giveFeed('Manifest saved. Total time elapsed:');
toc
