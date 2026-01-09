%% run_all_pupil_pe.m
%
% Batch script to generate pupil PE figures for all available sessions.
% Skips non-AV sessions and continues processing on errors.
%
% Author: Claude
% Date: 2026-01-09

function run_all_pupil_pe()

%% Setup Paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Define Session Paths
data_root = fullfile(findOneDrive, 'Neuronal Data Analysis');

% Find all session_data.mat files
session_files = dir(fullfile(data_root, '**', '*session_data.mat'));

fprintf('Found %d session files.\n', length(session_files));
fprintf('==========================================\n\n');

%% Process Each Session
n_success = 0;
n_skipped = 0;
n_failed = 0;
results_summary = {};

for i_sess = 1:length(session_files)
    file_path = fullfile(session_files(i_sess).folder, ...
        session_files(i_sess).name);

    % Extract session name
    [~, fname, ~] = fileparts(file_path);
    session_name = strrep(fname, '_session_data', '');

    fprintf('[%d/%d] Processing: %s\n', i_sess, length(session_files), ...
        session_name);

    try
        % Load session data
        loaded = load(file_path, 'session_data');
        session_data = loaded.session_data;

        % Get unique_id
        if isfield(session_data, 'metadata') && ...
                isfield(session_data.metadata, 'unique_id')
            unique_id = session_data.metadata.unique_id;
        else
            unique_id = session_name;
        end

        % Define conditions and check for AV session
        [conditions, is_av_session] = define_task_conditions( ...
            session_data.trialInfo, ...
            session_data.eventTimes, ...
            unique_id);

        if ~is_av_session
            fprintf('  -> SKIPPED: Not an AV session\n\n');
            n_skipped = n_skipped + 1;
            results_summary{end+1} = sprintf('%s: SKIPPED (not AV)', ...
                session_name);
            continue;
        end

        % Prepare core data
        if isfield(session_data, 'analysis') && ...
                isfield(session_data.analysis, 'core_data') && ...
                isfield(session_data.analysis.core_data, 'pupil') && ...
                isfield(session_data.analysis.core_data.pupil, 'outcomeOn')
            core_data = session_data.analysis.core_data;
        else
            if isfield(session_data, 'analysis') && ...
                    isfield(session_data.analysis, 'selected_neurons')
                selected_neurons = session_data.analysis.selected_neurons;
            else
                n_clusters = height(session_data.spikes.cluster_info);
                selected_neurons = false(n_clusters, 1);
            end
            core_data = prepare_core_data(session_data, ...
                selected_neurons, {'outcomeOn'});
        end

        % Run analysis
        pe_results = analyze_pupil_pe(core_data, conditions);

        % Generate plot
        fig = plot_pupil_pe(pe_results, unique_id);
        close(fig);

        fprintf('  -> SUCCESS\n\n');
        n_success = n_success + 1;
        results_summary{end+1} = sprintf('%s: SUCCESS', session_name);

    catch ME
        fprintf('  -> FAILED: %s\n\n', ME.message);
        n_failed = n_failed + 1;
        results_summary{end+1} = sprintf('%s: FAILED (%s)', ...
            session_name, ME.message);
    end
end

%% Print Summary
fprintf('\n==========================================\n');
fprintf('BATCH PROCESSING COMPLETE\n');
fprintf('==========================================\n');
fprintf('  Successful: %d\n', n_success);
fprintf('  Skipped:    %d (non-AV sessions)\n', n_skipped);
fprintf('  Failed:     %d\n', n_failed);
fprintf('==========================================\n\n');

fprintf('Session Results:\n');
for i = 1:length(results_summary)
    fprintf('  %s\n', results_summary{i});
end

end
