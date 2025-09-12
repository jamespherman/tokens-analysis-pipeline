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

function [aggregated_sc_data, aggregated_snc_data] = aggregate_analysis_results(manifest)

% start timer and define in-line function to give user feedback:
tic;
giveFeed = @(x)disp([num2str(round(toc, 1)) 's - ' x]);

%% Setup Paths
% Add the 'utils' directory to the path so that helper functions can be
% found.
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Discover Superset of Analysis Fields and Dimensions
% To decouple aggregation from a static analysis plan, we first iterate
% through all completed session files to dynamically discover the full
% superset of analysis results and key dimensions (e.g., n_time_bins).
all_roc_fields = {};
roc_time_bins = containers.Map('KeyType', 'char', 'ValueType', 'any');
discovered_anova_fields = struct();
n_time_bins_canonical = NaN; % To store the canonical number of time bins for ANOVA

complete_sessions = manifest(strcmp(manifest.analysis_status, 'complete'), :);

nSessions = height(complete_sessions);

% give user feedback:
giveFeed('Initial data loop to discover analysis fields & dimensions.')

% loop
for i_session = 1:nSessions

    % give user feedback:
    giveFeed(['Loading session ' num2str(i_session) ' of ' ...
        num2str(nSessions) '.'])

    % load data:
    session_id = complete_sessions.unique_id{i_session};
    session_data_path = fullfile(findOneDrive(), ...
        'Neuronal Data Analysis', session_id, ...
        [session_id '_session_data.mat']);

    % 
    giveFeed(['Session ' num2str(i_session) ' of ' ...
        num2str(nSessions) ' loaded.'])

    if exist(session_data_path, 'file')
        data = load(session_data_path, 'session_data');
        session_data = data.session_data;

        if isfield(session_data, 'analysis')
            % Discover ROC comparison fields and their time dimensions
            if isfield(session_data.analysis, 'roc_comparison')
                roc_fields = fieldnames(session_data.analysis.roc_comparison);
                all_roc_fields = [all_roc_fields; roc_fields];
                for i_field = 1:length(roc_fields)
                    field = roc_fields{i_field};
                    if ~isKey(roc_time_bins, field)
                        n_bins = size(session_data.analysis.roc_comparison.(field).sig, 2);
                        roc_time_bins(field) = n_bins;
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

                        % Capture canonical n_time_bins from the first valid p-value matrix
                        if isnan(n_time_bins_canonical) && ~isempty(p_value_fields)
                            p_field = p_value_fields{1};
                            n_time_bins_canonical = size(anova_results.(event_name).(p_field), 2);
                        end
                    end
                end
            end
        end
    end
end

% Get the unique field names for each analysis type
discovered_roc_fields = unique(all_roc_fields);

% Make the list of p-value fields for each ANOVA event unique
alignment_events = fieldnames(discovered_anova_fields);
for i_event = 1:length(alignment_events)
    event_name = alignment_events{i_event};
    discovered_anova_fields.(event_name) = unique(discovered_anova_fields.(event_name));
end


% The analysis plan is now discovered dynamically from the data,
% so the static define_analysis_plan() is no longer needed.

%% Initialize Aggregated Data Structures
aggregated_sc_data = struct();
aggregated_snc_data = struct();
brain_areas = {'SC', 'SNc'};

%% Main Loop

% give user feedback:
giveFeed('Main data loop over ''brain areas'' to aggregate.')

% loop
for i_area = 1:length(brain_areas)
    current_area = brain_areas{i_area};

    % Initialize a temporary struct for the current area's aggregated data
    aggregated_data = struct();

    % --- Dynamic Initialization from Discovered Fields ---
    % Initialize fields for ROC comparisons based on the discovered superset
    for i_comp = 1:length(discovered_roc_fields)
        comp_name = discovered_roc_fields{i_comp};
        aggregated_data.roc_comparison.(comp_name).sig = [];
    end

    % Initialize fields for ANOVA results using the discovered nested structure
    alignment_events = fieldnames(discovered_anova_fields);
    for i_event = 1:length(alignment_events)
        event_name = alignment_events{i_event};
        p_value_fields = discovered_anova_fields.(event_name);
        for i_p_field = 1:length(p_value_fields)
            p_field_name = p_value_fields{i_p_field};
            aggregated_data.anova_results.(event_name).(p_field_name) = [];
        end
    end


    % Filter manifest for complete sessions in the current area
    area_sessions = manifest(strcmp(manifest.brain_area, current_area) & ...
                             strcmp(manifest.analysis_status, 'complete'), :);

    % give user feedback:
    giveFeed(['Loop over ' brain_areas{i_area} ' sessions'])

    % Session Loop
    for i_session = 1:height(area_sessions)
        session_id = area_sessions.unique_id{i_session};

        % Construct path and load session data
        session_data_path = fullfile(findOneDrive(), ...
            'Neuronal Data Analysis', session_id, ...
            [session_id '_session_data.mat']);

        if ~exist(session_data_path, 'file')
            warning('aggregate_analysis_results:fileNotFound', ...
                    'Could not find session_data.mat for %s. Skipping.', ...
                    session_id);
            continue;
        end

        data = load(session_data_path, 'session_data');
        session_data = data.session_data;

        % Skip if analysis results are not present
        if ~isfield(session_data, 'analysis')
            warning('aggregate_analysis_results:noAnalysis', ...
                'No analysis field in session_data for %s. Skipping.', ...
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

        % give user feedback:
        giveFeed('Loop over discovered ROC fields.')

        % --- Aggregate ROC Comparison Results ---
        for i_roc = 1:length(discovered_roc_fields)
            field = discovered_roc_fields{i_roc};

            if isfield(session_data.analysis, 'roc_comparison') && ...
               isfield(session_data.analysis.roc_comparison, field)
                data_to_append = session_data.analysis.roc_comparison.(field).sig;
            else
                % Pad with NaNs of the correct size for the specific comparison
                n_bins = roc_time_bins(field);
                data_to_append = nan(n_neurons, n_bins);
            end

            aggregated_data.roc_comparison.(field).sig = ...
                [aggregated_data.roc_comparison.(field).sig; data_to_append];
        end

        % give user feedback:
        giveFeed('Loop over nested ANOVA results.')

        % --- Aggregate Nested ANOVA Results ---
        alignment_events = fieldnames(discovered_anova_fields);
        for i_event = 1:length(alignment_events)
            event_name = alignment_events{i_event};
            p_value_fields = discovered_anova_fields.(event_name);

            for i_p_field = 1:length(p_value_fields)
                p_field_name = p_value_fields{i_p_field};

                % Check if the specific p-value field exists for the session
                if isfield(session_data.analysis, 'anova_results') && ...
                   isfield(session_data.analysis.anova_results, event_name) && ...
                   isfield(session_data.analysis.anova_results.(event_name), p_field_name)
                    data_to_append = session_data.analysis.anova_results.(event_name).(p_field_name);
                else
                    % Pad with NaNs of the canonical size.
                    data_to_append = nan(n_neurons, n_time_bins_canonical);
                end
                try
                % Append data to the correct nested field
                aggregated_data.anova_results.(event_name).(p_field_name) = ...
                    [aggregated_data.anova_results.(event_name).(p_field_name); data_to_append];
                catch me
                    keyboard
                end
            end
        end
    end

    % Assign the populated temporary struct to the correct output variable
    if strcmp(current_area, 'SC')
        aggregated_sc_data = aggregated_data;
    else
        aggregated_snc_data = aggregated_data;
    end
end

end
