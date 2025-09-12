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

% Get the definitive analysis plan
analysis_plan = define_analysis_plan();

%% Initialize Aggregated Data Structures
aggregated_sc_data = struct();
aggregated_snc_data = struct();
brain_areas = {'SC', 'SNc'};

%% Main Loop
for i_area = 1:length(brain_areas)
    current_area = brain_areas{i_area};

    % Initialize a temporary struct for the current area's aggregated data
    aggregated_data = struct();

    % --- Dynamic Initialization from Analysis Plan ---
    % Initialize fields for ROC comparisons
    for i_comp = 1:length(analysis_plan.roc_comparison.comparisons_to_run)
        comp_name = analysis_plan.roc_comparison.comparisons_to_run(i_comp).name;
        aggregated_data.roc_comparison.(comp_name).sig = [];
    end

    % Initialize fields for ANOVA results
    if analysis_plan.anova.run
        for i_field = 1:length(analysis_plan.anova.fields_to_aggregate)
            field_name = analysis_plan.anova.fields_to_aggregate{i_field};
            aggregated_data.anova_results.(field_name) = [];
        end
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
        roc_comparisons = analysis_plan.roc_comparison.comparisons_to_run;
        for i_roc = 1:length(roc_comparisons)
            field = roc_comparisons(i_roc).name;

            % The field is guaranteed to exist in aggregated_data due to
            % the dynamic initialization step.
            aggregated_data.roc_comparison.(field).sig = [aggregated_data.roc_comparison.(field).sig; ...
                session_data.analysis.roc_comparison.(field).sig];
        end

        % --- Aggregate ANOVA Results (Flexible Handling) ---
        anova_results = session_data.analysis.anova_results;
        if analysis_plan.anova.run
            anova_fields = analysis_plan.anova.fields_to_aggregate;
            for i_anova = 1:length(anova_fields)
                field = anova_fields{i_anova};

                if isfield(anova_results, field)
                    data_to_append = anova_results.(field);
                else
                    % This handles cases where a field is not applicable to the
                    % session (e.g., flicker terms in a 'main' task session).
                    % We get nTimeBins from a field that is guaranteed to exist.
                    any_existing_field = fieldnames(anova_results);
                    any_existing_subfield = fieldnames(anova_results.(any_existing_field{1}));
                    n_time_bins = size(anova_results.(any_existing_field{1}).(any_existing_subfield{1}), 2);
                    data_to_append = nan(n_neurons, n_time_bins);
                end
                aggregated_data.anova_results.(field) = ...
                    [aggregated_data.anova_results.(field); data_to_append];
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
