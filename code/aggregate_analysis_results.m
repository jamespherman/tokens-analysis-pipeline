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
all_anova_fields = {};
n_time_bins_canonical = NaN; % To store the canonical number of time bins for ANOVA

complete_sessions = manifest(strcmp(manifest.analysis_status, 'complete'), :);

for i_session = 1:height(complete_sessions)
    session_id = complete_sessions.unique_id{i_session};
    session_data_path = fullfile(findOneDrive(), ...
        'Neuronal Data Analysis', session_id, ...
        [session_id '_session_data.mat']);

    if exist(session_data_path, 'file')
        data = load(session_data_path, 'session_data');
        session_data = data.session_data;

        if isfield(session_data, 'analysis')
            % Discover ROC comparison fields
            if isfield(session_data.analysis, 'roc_comparison')
                all_roc_fields = [all_roc_fields; fieldnames(session_data.analysis.roc_comparison)];
            end
            % Discover ANOVA fields and n_time_bins
            if isfield(session_data.analysis, 'anova_results')
                anova_results = session_data.analysis.anova_results;
                all_anova_fields = [all_anova_fields; fieldnames(anova_results)];

                % If we haven't found n_time_bins yet, get it from the first
                % available anova result in this session.
                if isnan(n_time_bins_canonical) && ~isempty(fieldnames(anova_results))
                    any_existing_field = fieldnames(anova_results);
                    any_existing_subfield = fieldnames(anova_results.(any_existing_field{1}));
                    n_time_bins_canonical = size(anova_results.(any_existing_field{1}).(any_existing_subfield{1}), 2);
                end
            end
        end
    end
end

% Get the unique field names for each analysis type
discovered_roc_fields = unique(all_roc_fields);
discovered_anova_fields = unique(all_anova_fields);


% The analysis plan is now discovered dynamically from the data,
% so the static define_analysis_plan() is no longer needed.

%% Initialize Aggregated Data Structures
aggregated_sc_data = struct();
aggregated_snc_data = struct();
brain_areas = {'SC', 'SNc'};

%% Main Loop
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

    % Initialize fields for ANOVA results based on the discovered superset
    for i_field = 1:length(discovered_anova_fields)
        field_name = discovered_anova_fields{i_field};
        aggregated_data.anova_results.(field_name) = [];
    end


    % Filter manifest for complete sessions in the current area
    area_sessions = manifest(strcmp(manifest.brain_area, current_area) & ...
                             strcmp(manifest.analysis_status, 'complete'), :);

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

        % --- Aggregate ROC Comparison Results ---
        for i_roc = 1:length(discovered_roc_fields)
            field = discovered_roc_fields{i_roc};

            if isfield(session_data.analysis, 'roc_comparison') && ...
               isfield(session_data.analysis.roc_comparison, field)
                data_to_append = session_data.analysis.roc_comparison.(field).sig;
            else
                % Pad with NaNs for missing data
                data_to_append = nan(n_neurons, 1);
            end

            aggregated_data.roc_comparison.(field).sig = ...
                [aggregated_data.roc_comparison.(field).sig; data_to_append];
        end

        % --- Aggregate ANOVA Results (Flexible Handling) ---
        for i_anova = 1:length(discovered_anova_fields)
            field = discovered_anova_fields{i_anova};

            if isfield(session_data.analysis, 'anova_results') && ...
               isfield(session_data.analysis.anova_results, field)
                data_to_append = session_data.analysis.anova_results.(field);
            else
                % Pad with NaNs of the canonical size.
                % If n_time_bins_canonical is still NaN (e.g., no ANOVA
                % results found anywhere), this will create a 0-column
                % matrix, which is handled but should be noted.
                data_to_append = nan(n_neurons, n_time_bins_canonical);
            end
            aggregated_data.anova_results.(field) = ...
                [aggregated_data.anova_results.(field); data_to_append];
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
