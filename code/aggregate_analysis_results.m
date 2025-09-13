%% aggregate_analysis_results.m
%
% Systematically gathers all pre-computed single-session analysis results
% and pools them by brain area. This creates the final data structures
% needed for population-level statistics and visualization.
%
% INPUT:
%   manifest - A table containing session metadata, including unique_id,
%              brain_area, and analysis_status.
%
% OUTPUT:
%   aggregated_sc_data  - A struct containing aggregated data for SC.
%   aggregated_snc_data - A struct containing aggregated data for SNc.
%
% Author: Jules
% Date: 2025-09-12
%

function [aggregated_sc_data, aggregated_snc_data] = aggregate_analysis_results

% start timer and define in-line function to give user feedback:
tic;
giveFeed = @(x)disp([num2str(round(toc, 1)) 's - ' x]);

%% Setup Paths
% Add the 'utils' directory to the path so that helper functions can be
% found.
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir); % Assumes script is in code/
addpath(project_root); % Add project root to path

%% Load Manifest
giveFeed('Loading session manifest...');
manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');
manifest = readtable(manifest_path);
giveFeed('Manifest loaded.');

%% Single Pass: Discover, Load, and Cache
% To improve efficiency, we iterate through all completed session files
% just once. In this single pass, we dynamically discover the superset of
% analysis fields, load the data from disk, and cache it in memory.
all_roc_fields = {};
roc_time_vectors = containers.Map('KeyType', 'char', 'ValueType', 'any');
all_baseline_comparison_fields = {};
baseline_comparison_time_vectors = containers.Map('KeyType', 'char', 'ValueType', 'any');
discovered_anova_fields = struct();
nBins = struct(); % Initialize a single, nested struct to hold all bin counts

complete_sessions = manifest(strcmp(manifest.analysis_status, 'complete'), :);
nSessions = height(complete_sessions);

% Initialize a cell array to cache the loaded session data
session_data_cache = cell(nSessions, 1);

giveFeed('Starting single pass to discover, load, and cache session data.');

% Loop through all completed sessions
for i_session = 1:nSessions
    session_id = complete_sessions.unique_id{i_session};
    giveFeed(['Processing session ' num2str(i_session) ' of ' ...
        num2str(nSessions) ': ' session_id]);

    % Construct path to session data
    session_data_path = fullfile(findOneDrive(), ...
        'Neuronal Data Analysis', session_id, ...
        [session_id '_session_data.mat']);

    if exist(session_data_path, 'file')
        % Load data from disk
        data = load(session_data_path, 'session_data');
        session_data = data.session_data;

        % --- Cache the loaded data ---
        session_data_cache{i_session} = session_data;

        % --- "Discovery" Logic: Dynamically find all analysis fields ---
        if isfield(session_data, 'analysis')
            % Discover ROC comparison fields and their time dimensions
            if isfield(session_data.analysis, 'roc_comparison')
                roc_fields = fieldnames(session_data.analysis.roc_comparison);
                all_roc_fields = [all_roc_fields; roc_fields];
                for i_field = 1:length(roc_fields)
                    field = roc_fields{i_field};
                    if ~isfield(nBins, 'roc_comparison') || ~isfield(nBins.roc_comparison, field)
                        n_bins = size(session_data.analysis.roc_comparison.(field).sig, 2);
                        nBins.roc_comparison.(field) = n_bins;
                        % Also store the time vector for this comparison
                        roc_time_vectors(field) = session_data.analysis.roc_comparison.(field).time_vector;
                    end
                end
            end

            % Discover baseline_comparison fields and their time dimensions
            if isfield(session_data.analysis, 'baseline_comparison')
                baseline_comp_fields = fieldnames(session_data.analysis.baseline_comparison);
                all_baseline_comparison_fields = [all_baseline_comparison_fields; baseline_comp_fields];
                for i_field = 1:length(baseline_comp_fields)
                    field = baseline_comp_fields{i_field};
                    if ~isfield(nBins, 'baseline_comparison') || ~isfield(nBins.baseline_comparison, field)
                        n_bins = size(session_data.analysis.baseline_comparison.(field).sig, 2);
                        nBins.baseline_comparison.(field) = n_bins;
                        % Also store the time vector for this comparison
                        baseline_comparison_time_vectors(field) = session_data.analysis.baseline_comparison.(field).time_vector;
                    end
                end
            end

            % Discover nested ANOVA fields (alignment_event -> p_value_field)
            if isfield(session_data.analysis, 'anova_results')
                anova_results = session_data.analysis.anova_results;
                alignment_events = fieldnames(anova_results);

                for i_event = 1:length(alignment_events)
                    event_name = alignment_events{i_event};
                    if isstruct(anova_results.(event_name))
                        p_value_fields = fieldnames(anova_results.(event_name));

                        if ~isfield(discovered_anova_fields, event_name)
                            discovered_anova_fields.(event_name) = {};
                        end
                        discovered_anova_fields.(event_name) = [discovered_anova_fields.(event_name); p_value_fields];

                        % Store n_time_bins for each analysis
                        for i_p_field = 1:length(p_value_fields)
                            p_field = p_value_fields{i_p_field};
                            if ~isfield(nBins, 'anova_results') || ...
                               ~isfield(nBins.anova_results, event_name) || ...
                               ~isfield(nBins.anova_results.(event_name), p_field)
                                n_bins = size(anova_results.(event_name).(p_field), 2);
                                nBins.anova_results.(event_name).(p_field) = n_bins;
                            end
                        end
                    end
                end
            end
        end
    else
        warning('aggregate_analysis_results:fileNotFound', ...
                'Could not find session_data.mat for %s. Skipping.', ...
                session_id);
        % Store an empty value in the cache to maintain indexing
        session_data_cache{i_session} = [];
    end
end

% Get the unique field names for each analysis type
discovered_roc_fields = unique(all_roc_fields);
discovered_baseline_comparison_fields = unique(all_baseline_comparison_fields);

% Make the list of p-value fields for each ANOVA event unique
alignment_events = fieldnames(discovered_anova_fields);
for i_event = 1:length(alignment_events)
    event_name = alignment_events{i_event};
    discovered_anova_fields.(event_name) = unique(discovered_anova_fields.(event_name));
end

%% Initialize Aggregated Data Structures
aggregated_sc_data = struct();
aggregated_snc_data = struct();
brain_areas = {'SC', 'SNc'};

%% Aggregation Pass (from cached data)
giveFeed('Aggregating data from in-memory cache.');

for i_area = 1:length(brain_areas)
    current_area = brain_areas{i_area};
    aggregated_data = struct();

    % --- Dynamic Initialization from Discovered Fields ---
    for i_comp = 1:length(discovered_roc_fields)
        comp_name = discovered_roc_fields{i_comp};
        aggregated_data.roc_comparison.(comp_name).sig = [];
        aggregated_data.roc_comparison.(comp_name).time_vector = roc_time_vectors(comp_name);
    end

    for i_comp = 1:length(discovered_baseline_comparison_fields)
        comp_name = discovered_baseline_comparison_fields{i_comp};
        aggregated_data.baseline_comparison.(comp_name).sig = [];
        aggregated_data.baseline_comparison.(comp_name).time_vector = baseline_comparison_time_vectors(comp_name);
    end

    alignment_events = fieldnames(discovered_anova_fields);
    for i_event = 1:length(alignment_events)
        event_name = alignment_events{i_event};
        p_value_fields = discovered_anova_fields.(event_name);
        for i_p_field = 1:length(p_value_fields)
            p_field_name = p_value_fields{i_p_field};
            aggregated_data.anova_results.(event_name).(p_field_name) = [];
        end
    end

    % Find indices of sessions corresponding to the current brain area
    area_session_indices = find(strcmp(complete_sessions.brain_area, current_area));

    giveFeed(['Aggregating ' num2str(length(area_session_indices)) ...
        ' sessions for ' current_area]);

    % Loop through the sessions for the current brain area
    for i_idx = 1:length(area_session_indices)
        session_idx_in_cache = area_session_indices(i_idx);

        % --- Retrieve session data from the cache ---
        session_data = session_data_cache{session_idx_in_cache};

        % Skip if data was not loaded or is empty
        if isempty(session_data)
            continue;
        end

        session_id = complete_sessions.unique_id{session_idx_in_cache};

        if ~isfield(session_data, 'analysis')
            warning('aggregate_analysis_results:noAnalysis', ...
                'No analysis field in cached data for %s. Skipping.', ...
                session_id);
            continue;
        end

        n_neurons = size(session_data.analysis.selected_neurons, 1);

        % --- Aggregate Session ID ---
        session_ids_to_append = repmat({session_id}, n_neurons, 1);
        if ~isfield(aggregated_data, 'session_id')
            aggregated_data.session_id = session_ids_to_append;
        else
            aggregated_data.session_id = [aggregated_data.session_id; session_ids_to_append];
        end

        % --- Aggregate ROC Comparison Results ---
        for i_roc = 1:length(discovered_roc_fields)
            field = discovered_roc_fields{i_roc};
            if isfield(session_data.analysis.roc_comparison, field)
                data_to_append = session_data.analysis.roc_comparison.(field).sig;
            else
                n_bins = nBins.roc_comparison.(field);
                data_to_append = nan(n_neurons, n_bins);
            end
            aggregated_data.roc_comparison.(field).sig = ...
                [aggregated_data.roc_comparison.(field).sig; data_to_append];
        end

        % --- Aggregate Baseline Comparison Results ---
        for i_bc = 1:length(discovered_baseline_comparison_fields)
            field = discovered_baseline_comparison_fields{i_bc};
            if isfield(session_data.analysis, 'baseline_comparison') && isfield(session_data.analysis.baseline_comparison, field)
                data_to_append = session_data.analysis.baseline_comparison.(field).sig;
            else
                n_bins = nBins.baseline_comparison.(field);
                data_to_append = nan(n_neurons, n_bins);
            end
            aggregated_data.baseline_comparison.(field).sig = ...
                [aggregated_data.baseline_comparison.(field).sig; data_to_append];
        end

        % --- Aggregate Nested ANOVA Results ---
        for i_event = 1:length(alignment_events)
            event_name = alignment_events{i_event};
            p_value_fields = discovered_anova_fields.(event_name);
            for i_p_field = 1:length(p_value_fields)
                p_field_name = p_value_fields{i_p_field};
                if isfield(session_data.analysis.anova_results, event_name) && ...
                   isfield(session_data.analysis.anova_results.(event_name), p_field_name)
                    data_to_append = session_data.analysis.anova_results.(event_name).(p_field_name);
                else
                    n_bins = nBins.anova_results.(event_name).(p_field_name);
                    data_to_append = nan(n_neurons, n_bins);
                end
                aggregated_data.anova_results.(event_name).(p_field_name) = ...
                    [aggregated_data.anova_results.(event_name).(p_field_name); data_to_append];
            end
        end
    end

    % Assign the populated struct to the correct output variable
    if strcmp(current_area, 'SC')
        aggregated_sc_data = aggregated_data;
    else
        aggregated_snc_data = aggregated_data;
    end
end

% save the data:
giveFeed('Saving aggregated data to file...');
saveFileName = fullfile(project_root, 'data', 'processed', ...
    'aggregated_analysis_data.mat');
save(saveFileName, 'aggregated_sc_data', 'aggregated_snc_data');
giveFeed('Aggregation complete.');

end
