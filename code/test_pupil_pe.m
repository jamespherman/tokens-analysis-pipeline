%% test_pupil_pe.m
%
% Test script for verifying the pupil PE (RPE/SPE/Interaction) pipeline on
% a single session. Loads session data, checks for AV trials, prepares
% core data, runs the analysis, generates plots, and prints diagnostics.
%
% USAGE:
%   test_pupil_pe()          - Prompts user to select file via uigetfile
%   test_pupil_pe(file_path) - Uses the provided path to session_data.mat
%
% The script will error if the selected session does not contain AV trials,
% as the pupil PE analysis requires SPE conditions.
%
% Author: Claude
% Date: 2026-01-09

function test_pupil_pe(file_path)

%% Setup Paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

tic;
giveFeed = @(x)disp([num2str(toc) 's - ' x]);

%% Step 1: Load Session Data
giveFeed('Step 1: Loading session data...');

if nargin < 1 || isempty(file_path)
    [filename, pathname] = uigetfile('*.mat', ...
        'Select a session_data.mat file');

    if isequal(filename, 0)
        error('test_pupil_pe:userCancelled', ...
            'File selection cancelled by user.');
    end

    file_path = fullfile(pathname, filename);
end

if ~exist(file_path, 'file')
    error('test_pupil_pe:fileNotFound', ...
        'File not found: %s', file_path);
end

fprintf('Loading: %s\n', file_path);
loaded_data = load(file_path, 'session_data');

if ~isfield(loaded_data, 'session_data')
    error('test_pupil_pe:invalidFile', ...
        'File does not contain a ''session_data'' variable.');
end

session_data = loaded_data.session_data;

if isfield(session_data, 'metadata') && ...
        isfield(session_data.metadata, 'unique_id')
    unique_id = session_data.metadata.unique_id;
else
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
    error('test_pupil_pe:conditionError', ...
        'Failed to define task conditions: %s', ME.message);
end

if ~is_av_session
    error('test_pupil_pe:notAVSession', ...
        ['This session has no AV trials. Select an AV session.\n' ...
        'The pupil PE analysis requires SPE conditions.']);
end

giveFeed('Session confirmed as AV session.');

%% Step 3: Prepare Core Data
giveFeed('Step 3: Preparing core data...');

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

if ~has_valid_core_data
    if isfield(session_data, 'analysis') && ...
            isfield(session_data.analysis, 'selected_neurons')
        selected_neurons = session_data.analysis.selected_neurons;
    else
        n_clusters = height(session_data.spikes.cluster_info);
        selected_neurons = false(n_clusters, 1);
        giveFeed('No selected_neurons found. Using empty selection.');
    end

    alignment_events = {'outcomeOn'};

    try
        core_data = prepare_core_data(session_data, selected_neurons, ...
            alignment_events);
    catch ME
        error('test_pupil_pe:prepareDataError', ...
            'Failed to prepare core data: %s', ME.message);
    end
end

if ~isfield(core_data, 'pupil') || ...
        ~isfield(core_data.pupil, 'outcomeOn')
    error('test_pupil_pe:missingPupilData', ...
        'core_data does not contain pupil.outcomeOn data.');
end

giveFeed(sprintf('Core data prepared. Pupil traces: %d trials x %d samples', ...
    size(core_data.pupil.outcomeOn.traces, 1), ...
    size(core_data.pupil.outcomeOn.traces, 2)));

%% Step 4: Run Pupil PE Analysis
giveFeed('Step 4: Running pupil PE analysis...');

try
    results = analyze_pupil_pe(core_data, conditions);
catch ME
    error('test_pupil_pe:analysisError', ...
        'Failed to run pupil PE analysis: %s', ME.message);
end

giveFeed('Pupil PE analysis complete.');

%% Step 5: Generate Plot
giveFeed('Step 5: Generating pupil PE plot...');

try
    fig = plot_pupil_pe(results, unique_id);
catch ME
    error('test_pupil_pe:plotError', ...
        'Failed to generate pupil PE plot: %s', ME.message);
end

giveFeed('Plot generated successfully.');

%% Step 6: Print Diagnostic Output
giveFeed('Step 6: Printing diagnostic summary...');

fprintf('\n');
fprintf('============================================================\n');
fprintf('              PUPIL PE DIAGNOSTIC SUMMARY                   \n');
fprintf('============================================================\n');
fprintf('Session: %s\n', unique_id);
fprintf('\n');

alpha = 0.05;

% Panel A: RPE
fprintf('--- Panel A: RPE (Normal Distribution) ---\n');
fprintf('  Trial counts: Rare-Low=%d, Common=%d, Rare-High=%d\n', ...
    results.rpe.trial_counts(1), results.rpe.trial_counts(2), ...
    results.rpe.trial_counts(3));

n_sig_rpe = sum(results.rpe.p_rpe < alpha, 'omitnan');
n_total = sum(~isnan(results.rpe.p_rpe));
fprintf('  Significant bins (p < %.2f): %d / %d (%.1f%%)\n', ...
    alpha, n_sig_rpe, n_total, 100 * n_sig_rpe / n_total);

[min_p, idx] = min(results.rpe.p_rpe);
if ~isnan(min_p)
    fprintf('  Peak effect: t = %+.3f s, p = %.4f\n', ...
        results.time_vector(idx), min_p);
end
fprintf('\n');

% Panel B: SPE
fprintf('--- Panel B: SPE ---\n');
fprintf('  Trial counts: No-Flicker=%d, Omitted=%d, Surprising=%d, ', ...
    results.spe.trial_counts(1), results.spe.trial_counts(2), ...
    results.spe.trial_counts(3));
fprintf('Certain=%d\n', results.spe.trial_counts(4));

n_sig_spe = sum(results.spe.p_contrast < alpha, 'omitnan');
n_total = sum(~isnan(results.spe.p_contrast));
fprintf('  Significant contrast bins (p < %.2f): %d / %d (%.1f%%)\n', ...
    alpha, n_sig_spe, n_total, 100 * n_sig_spe / n_total);

[min_p, idx] = min(results.spe.p_contrast);
if ~isnan(min_p)
    fprintf('  Peak contrast: t = %+.3f s, p = %.4f\n', ...
        results.time_vector(idx), min_p);
end
fprintf('\n');

% Panel C: Interaction
fprintf('--- Panel C: RPE x AV Interaction ---\n');
fprintf('  Trial counts:\n');
fprintf('                No-AV    AV\n');
fprintf('    Common:     %4d   %4d\n', ...
    results.interaction.trial_counts(1, 1), ...
    results.interaction.trial_counts(1, 2));
fprintf('    Rare-High:  %4d   %4d\n', ...
    results.interaction.trial_counts(2, 1), ...
    results.interaction.trial_counts(2, 2));

n_sig_rpe_main = sum(results.interaction.p_rpe < alpha, 'omitnan');
n_sig_av_main = sum(results.interaction.p_av < alpha, 'omitnan');
n_sig_int = sum(results.interaction.p_int < alpha, 'omitnan');
n_total = sum(~isnan(results.interaction.p_int));

fprintf('  Significant bins (p < %.2f):\n', alpha);
fprintf('    RPE main:      %3d / %d (%.1f%%)\n', ...
    n_sig_rpe_main, n_total, 100 * n_sig_rpe_main / n_total);
fprintf('    AV main:       %3d / %d (%.1f%%)\n', ...
    n_sig_av_main, n_total, 100 * n_sig_av_main / n_total);
fprintf('    Interaction:   %3d / %d (%.1f%%)\n', ...
    n_sig_int, n_total, 100 * n_sig_int / n_total);

[min_p, idx] = min(results.interaction.p_int);
if ~isnan(min_p)
    fprintf('  Peak interaction: t = %+.3f s, p = %.4f\n', ...
        results.time_vector(idx), min_p);
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
