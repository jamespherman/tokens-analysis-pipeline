%% test_pupil_anova.m
%
% Test script for verifying the pupil ANOVA pipeline on a single session.
% Loads session data, checks for AV trials, prepares core data, runs the
% ANOVA analysis, generates plots, and prints diagnostic output.
%
% USAGE:
%   test_pupil_anova()          - Prompts user to select file via uigetfile
%   test_pupil_anova(file_path) - Uses the provided path to session_data.mat
%
% The script will error if the selected session does not contain AV trials,
% as the pupil ANOVA analysis requires the AV probability factor.
%
% Author: Claude
% Date: 2026-01-09

function test_pupil_anova(file_path)

%% Setup Paths
% Add the 'utils' directory to the path so that helper functions can be
% found.
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

% Start timer
tic;

% In-line function to report timing
giveFeed = @(x)disp([num2str(toc) 's - ' x]);

%% Step 1: Load Session Data
giveFeed('Step 1: Loading session data...');

if nargin < 1 || isempty(file_path)
    % Prompt user to select a file
    [filename, pathname] = uigetfile('*.mat', ...
        'Select a session_data.mat file');

    if isequal(filename, 0)
        error('test_pupil_anova:userCancelled', ...
            'File selection cancelled by user.');
    end

    file_path = fullfile(pathname, filename);
end

% Check if file exists
if ~exist(file_path, 'file')
    error('test_pupil_anova:fileNotFound', ...
        'File not found: %s', file_path);
end

fprintf('Loading: %s\n', file_path);
loaded_data = load(file_path, 'session_data');

if ~isfield(loaded_data, 'session_data')
    error('test_pupil_anova:invalidFile', ...
        'File does not contain a ''session_data'' variable.');
end

session_data = loaded_data.session_data;

% Extract unique_id for titles and output
if isfield(session_data, 'metadata') && ...
        isfield(session_data.metadata, 'unique_id')
    unique_id = session_data.metadata.unique_id;
else
    % Fall back to filename
    [~, unique_id, ~] = fileparts(file_path);
    unique_id = strrep(unique_id, '_session_data', '');
end

giveFeed(sprintf('Loaded session: %s', unique_id));

%% Step 2: Define Task Conditions and Check for AV Session
giveFeed('Step 2: Defining task conditions...');

try
    [conditions, is_av_session] = define_task_conditions( ...
        session_data.trialInfo, ...
        session_data.eventTimes, ...
        unique_id);
catch ME
    error('test_pupil_anova:conditionError', ...
        'Failed to define task conditions: %s', ME.message);
end

% Check if this is an AV session
if ~is_av_session
    error('test_pupil_anova:notAVSession', ...
        ['This session has no AV trials. Select an AV session.\n' ...
        'The pupil ANOVA analysis requires AV probability as a factor.']);
end

giveFeed('Session confirmed as AV session.');

%% Step 3: Prepare Core Data
giveFeed('Step 3: Preparing core data...');

% Check if core_data is already prepared with outcomeOn alignment
has_valid_core_data = false;
if isfield(session_data, 'analysis') && ...
        isfield(session_data.analysis, 'core_data')
    core_data = session_data.analysis.core_data;
    if isfield(core_data, 'pupil') && ...
            isfield(core_data.pupil, 'outcomeOn')
        giveFeed('Using pre-computed core_data from session_data.analysis.');
        has_valid_core_data = true;
    else
        giveFeed(['Pre-computed core_data missing outcomeOn alignment. ' ...
            'Re-preparing...']);
    end
end

% Prepare core_data if needed
if ~has_valid_core_data
    % We need selected_neurons for prepare_core_data, but for pupil analysis
    % we can use an empty or minimal selection
    if isfield(session_data, 'analysis') && ...
            isfield(session_data.analysis, 'selected_neurons')
        selected_neurons = session_data.analysis.selected_neurons;
    else
        % For pupil-only analysis, create a dummy selection
        n_clusters = height(session_data.spikes.cluster_info);
        selected_neurons = false(n_clusters, 1);
        giveFeed('No selected_neurons found. Using empty selection.');
    end

    % Define alignment events - we need outcomeOn for pupil ANOVA
    alignment_events = {'outcomeOn'};

    try
        core_data = prepare_core_data(session_data, selected_neurons, ...
            alignment_events);
    catch ME
        error('test_pupil_anova:prepareDataError', ...
            'Failed to prepare core data: %s', ME.message);
    end
end

% Final verification that outcomeOn pupil data exists
if ~isfield(core_data, 'pupil') || ...
        ~isfield(core_data.pupil, 'outcomeOn')
    error('test_pupil_anova:missingPupilData', ...
        'core_data does not contain pupil.outcomeOn data.');
end

giveFeed(sprintf('Core data prepared. Pupil traces: %d trials x %d samples', ...
    size(core_data.pupil.outcomeOn.traces, 1), ...
    size(core_data.pupil.outcomeOn.traces, 2)));

%% Step 4: Run Pupil ANOVA Analysis
giveFeed('Step 4: Running pupil ANOVA analysis...');

try
    results = analyze_pupil_anova(core_data, conditions);
catch ME
    error('test_pupil_anova:analysisError', ...
        'Failed to run pupil ANOVA analysis: %s', ME.message);
end

giveFeed('Pupil ANOVA analysis complete.');

%% Step 5: Generate Plot
giveFeed('Step 5: Generating pupil ANOVA plot...');

try
    fig = plot_pupil_anova(results, unique_id);
catch ME
    error('test_pupil_anova:plotError', ...
        'Failed to generate pupil ANOVA plot: %s', ME.message);
end

giveFeed('Plot generated successfully.');

%% Step 6: Print Diagnostic Output
giveFeed('Step 6: Printing diagnostic summary...');

fprintf('\n');
fprintf('============================================================\n');
fprintf('              PUPIL ANOVA DIAGNOSTIC SUMMARY                \n');
fprintf('============================================================\n');
fprintf('Session: %s\n', unique_id);
fprintf('\n');

% Trial counts per cell (already printed by analyze_pupil_anova, but
% re-print for clarity)
fprintf('--- Trial Counts Per Cell ---\n');
fprintf('              AV 0%%    AV 50%%   AV 100%%\n');
fprintf('    Normal:   %4d      %4d      %4d\n', ...
    results.trial_counts(1, 1), results.trial_counts(1, 2), ...
    results.trial_counts(1, 3));
fprintf('    Uniform:  %4d      %4d      %4d\n', ...
    results.trial_counts(2, 1), results.trial_counts(2, 2), ...
    results.trial_counts(2, 3));
fprintf('\n');

% Count significant bins for each factor
alpha = 0.05;
n_sig_dist = sum(results.p_dist < alpha, 'omitnan');
n_sig_av = sum(results.p_av_prob < alpha, 'omitnan');
n_sig_int = sum(results.p_interaction < alpha, 'omitnan');
n_total_bins = sum(~isnan(results.p_dist));

fprintf('--- Significant Bins (p < %.2f) ---\n', alpha);
fprintf('  Distribution:    %3d / %d bins (%.1f%%)\n', ...
    n_sig_dist, n_total_bins, 100 * n_sig_dist / n_total_bins);
fprintf('  AV Probability:  %3d / %d bins (%.1f%%)\n', ...
    n_sig_av, n_total_bins, 100 * n_sig_av / n_total_bins);
fprintf('  Interaction:     %3d / %d bins (%.1f%%)\n', ...
    n_sig_int, n_total_bins, 100 * n_sig_int / n_total_bins);
fprintf('\n');

% Time of peak effect (lowest p-value) for each factor
fprintf('--- Peak Effect (Lowest p-value) ---\n');

% Distribution
[min_p_dist, idx_dist] = min(results.p_dist);
if ~isnan(min_p_dist)
    peak_time_dist = results.time_vector(idx_dist);
    fprintf('  Distribution:    t = %+.3f s, p = %.4f\n', ...
        peak_time_dist, min_p_dist);
else
    fprintf('  Distribution:    No valid p-values\n');
end

% AV Probability
[min_p_av, idx_av] = min(results.p_av_prob);
if ~isnan(min_p_av)
    peak_time_av = results.time_vector(idx_av);
    fprintf('  AV Probability:  t = %+.3f s, p = %.4f\n', ...
        peak_time_av, min_p_av);
else
    fprintf('  AV Probability:  No valid p-values\n');
end

% Interaction
[min_p_int, idx_int] = min(results.p_interaction);
if ~isnan(min_p_int)
    peak_time_int = results.time_vector(idx_int);
    fprintf('  Interaction:     t = %+.3f s, p = %.4f\n', ...
        peak_time_int, min_p_int);
else
    fprintf('  Interaction:     No valid p-values\n');
end

fprintf('\n');
fprintf('============================================================\n');
fprintf('Time window: [%.2f, %.2f] s relative to outcome onset\n', ...
    results.time_vector(1), results.time_vector(end));
fprintf('Bin width: 100 ms\n');
fprintf('============================================================\n');
fprintf('\n');

giveFeed('Test script complete.');

end
