%% aggregate_analysis_results.m
%
% Refactored to use a three-phase approach:
% 1. Single-pass load, discovery of dimensions, and caching.
% 2. Aggregation from the in-memory cache, with correct NaN-padding.
% 3. Saving the final aggregated data.
%
% This approach is robust to missing or heterogeneous analysis results
% across sessions.
%
% Author: Jules
% Date: 2025-09-14
%

function [aggregated_sc_data, aggregated_snc_data] = ...
    aggregate_analysis_results()

% Start timer and define in-line function to give user feedback:
tic;
giveFeed = @(x)disp([num2str(round(toc, 1)) 's - ' x]);

%% Phase 1: Initialization
giveFeed('Phase 1: Initializing...');

% Setup Paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);
addpath(project_root);

% Load Manifest
manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');
manifest = readtable(manifest_path);

% Initialize Caching and Discovery Structures
session_data_cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
dim_info = struct();

% Initialize Final Aggregated Data Structures
brain_areas = {'SC', 'SNc'};
aggregated_data_by_area = ...
    containers.Map('KeyType', 'char', 'ValueType', 'any');

for i_area = 1:length(brain_areas)
    area = brain_areas{i_area};
    agg_data = struct();
    agg_data.session_id = {};
    % Note: Analysis-specific fields are NOT initialized here.
    % They will be created dynamically during the aggregation phase based
    % on the discovered dimensions.
    aggregated_data_by_area(area) = agg_data;
end

giveFeed('Initialization complete.');

%% Phase 2: Single-Pass Data Loading, Discovery, and Caching
giveFeed('Phase 2: Loading, discovering, and caching session data...');
complete_sessions = manifest(strcmp(manifest.analysis_status, 'complete'), :);
nSessions = height(complete_sessions);
giveFeed(sprintf('Found %d completed sessions to process.', nSessions));

for i_session = 1:nSessions
    session_id = complete_sessions.unique_id{i_session};
    giveFeed(sprintf('Processing session %d/%d: %s', i_session, nSessions, session_id));

    session_data_path = fullfile(findOneDrive(), 'Neuronal Data Analysis', ...
        session_id, [session_id '_session_data.mat']);

    if ~exist(session_data_path, 'file')
        warning('aggregate_analysis_results:fileNotFound', ...
                'Could not find session_data.mat for %s. Skipping.', session_id);
        continue;
    end

    % Load the data
    data = load(session_data_path, 'session_data');
    session_data = data.session_data;

    % Cache the entire session_data struct
    session_data_cache(session_id) = session_data;

    % Discover dimensions if analysis data exists
    if isfield(session_data, 'analysis')
        dim_info = discover_dimensions_recursive(session_data.analysis, dim_info);
    else
        warning('aggregate_analysis_results:noAnalysis', ...
            'No analysis field in data for %s. Skipping discovery.', session_id);
    end
end
giveFeed('Data loading, discovery, and caching complete.');

%% Phase 3: Aggregation from Cache
giveFeed('Phase 3: Aggregating data from cache...');

for i_session = 1:nSessions
    session_id = complete_sessions.unique_id{i_session};
    brain_area = complete_sessions.brain_area{i_session};
    giveFeed(sprintf('Aggregating session %d/%d: %s', ...
        i_session, nSessions, session_id));

    % Retrieve data from cache
    session_data = session_data_cache(session_id);

    % Get the correct aggregated struct for the current brain area
    agg_data = aggregated_data_by_area(brain_area);

    % Determine the number of neurons for NaN padding
    n_neurons = 0;
    if isfield(session_data, 'analysis') && ...
       isfield(session_data.analysis, 'selected_neurons')
        n_neurons = size(session_data.analysis.selected_neurons, 1);
    end

    if n_neurons == 0
        giveFeed(sprintf('Skipping session %s: no selected neurons.', ...
            session_id));
        continue;
    end

    % Append session ID
    agg_data.session_id = [agg_data.session_id; ...
        repmat({session_id}, n_neurons, 1)];

    % Recursively aggregate all fields based on the dim_info template
    session_analysis_data = struct();
    if isfield(session_data, 'analysis')
        session_analysis_data = session_data.analysis;
    end

    agg_data = aggregate_recursive(dim_info, session_analysis_data, ...
        agg_data, n_neurons);

    % Update the map with the modified struct
    aggregated_data_by_area(brain_area) = agg_data;
end
giveFeed('Aggregation from cache complete.');

%% Finalize and Save
giveFeed('Saving aggregated data to file...');
aggregated_sc_data = aggregated_data_by_area('SC');
aggregated_snc_data = aggregated_data_by_area('SNc');

saveFileName = fullfile(project_root, 'data', 'processed', ...
    'aggregated_analysis_data.mat');
save(saveFileName, 'aggregated_sc_data', 'aggregated_snc_data', '-v7.3');
giveFeed('Aggregation complete.');

end

% =========================================================================
% Helper function to recursively discover dimensions
% =========================================================================
function dim_info = discover_dimensions_recursive(current_data, dim_info)
    fields = fieldnames(current_data);
    for i = 1:length(fields)
        field_name = fields{i};

        % Skip metadata fields that are not analysis results
        if ismember(field_name, ...
           {'selected_neurons', 'scSide', 'sig_epoch_comparison'})
            continue;
        end

        if isstruct(current_data.(field_name))
            if ~isfield(dim_info, field_name)
                dim_info.(field_name) = struct();
            end
            dim_info.(field_name) = discover_dimensions_recursive(...
                current_data.(field_name), dim_info.(field_name));
        else
            % This is a leaf node (a data matrix). Record its size.
            dims = size(current_data.(field_name));
            if ~isfield(dim_info, field_name)
                 % Store dimensions beyond the first (neuron) dimension
                 dim_info.(field_name).size = dims(2:end);
                 % Also store time vectors if they exist
                 if isfield(current_data, 'time_vector')
                     dim_info.(field_name).time_vector = ...
                         current_data.time_vector;
                 end
            end
        end
    end
end

% =========================================================================
% Helper function to recursively aggregate data
% =========================================================================
function agg_data = aggregate_recursive(dim_info_sub, ...
    session_analysis_sub, agg_data, n_neurons)

    dim_fields = fieldnames(dim_info_sub);
    for i = 1:length(dim_fields)
        field_name = dim_fields{i};
        current_dim_info = dim_info_sub.(field_name);

        % Check if this is a leaf node (an analysis result matrix)
        if isfield(current_dim_info, 'size')
            if ~isfield(agg_data, field_name)
                agg_data.(field_name) = [];
            end

            % Check if data exists in the current session
            if isfield(session_analysis_sub, field_name) && ...
               ~isempty(session_analysis_sub.(field_name))
                data_to_append = session_analysis_sub.(field_name);

                if isfield(current_dim_info, 'time_vector')
                    if ~isfield(agg_data, 'time_vectors') || ...
                       ~isfield(agg_data.time_vectors, field_name)
                        agg_data.time_vectors.(field_name) = ...
                            current_dim_info.time_vector;
                    end
                end
            else
                % Data is missing, so create a NaN matrix for padding
                pad_size = current_dim_info.size;
                if isempty(pad_size)
                    pad_size = [1]; % Handle scalar values
                end
                data_to_append = nan([n_neurons, pad_size]);
            end
            agg_data.(field_name) = [agg_data.(field_name); data_to_append];

        else % This is a struct, so we recurse deeper
            if ~isfield(agg_data, field_name)
                agg_data.(field_name) = struct();
            end

            next_session_sub = struct();
            if isfield(session_analysis_sub, field_name)
                next_session_sub = session_analysis_sub.(field_name);
            end

            agg_data.(field_name) = aggregate_recursive(...
                current_dim_info, ...
                next_session_sub, ...
                agg_data.(field_name), ...
                n_neurons);
        end
    end
end
