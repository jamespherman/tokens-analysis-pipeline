%% test_core_data_preparation.m
%
%   This script tests the refactored `prepare_core_data` module. It loads
%   a sample session, runs the data preparation pipeline, and generates a
%   figure with diagnostic plots to verify the output for all three key
%   alignment events: CUE_ON, outcomeOn, and reward.
%
%   The figure contains two rows:
%   - Top Row: Grand average PSTHs for neuronal data.
%   - Bottom Row: Grand average processed pupil traces.
%
% Author: Jules
% Date: 2025-09-08

%% Setup
clear; clc; close all;

% Start timer
tic;

% In-line function to report timing
giveFeed = @(x)disp([num2str(toc) 's - ' x]);

% Add utility functions to the MATLAB path
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Setup
unique_id = 'Feynman_08_07_2025_SNc'; % Hardcoded session
giveFeed(sprintf('Testing diagnostic workflow for session: %s', ...
    unique_id));

%% Load Manifest
giveFeed('Step 1: Loading session manifest...');
project_root = fullfile(findOneDrive, 'Code', 'tokens-analysis-pipeline');
manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');
manifest = readtable(manifest_path);

% Find the row for the current session
session_idx = find(strcmp(manifest.unique_id, unique_id));
if isempty(session_idx)
    error('test_core_data_preparation:sessionNotFound', ...
          'Session with unique_id "%s" not found in the manifest.', ...
          unique_id);
end

%% Load Data
% Define the project root and construct the path to a sample data file.
% NOTE: The user should replace 'sample_session.mat' with a real session
% file from their 'data/raw' directory.
giveFeed('Step 2: Loading session data...');
one_drive_path = findOneDrive;
data_file = fullfile(one_drive_path, ...
    'Neuronal Data Analysis', unique_id, ...
    [unique_id '_session_data.mat']);

% Check if the sample file exists
if ~exist(data_file, 'file')
    error('test_core_data_preparation:sampleNotFound', ...
          ['Sample data file not found. Please specify a valid path ' ...
          'to a \n session_data.mat file in the "Load Data" section of ' ...
          'this script.']);
end

fprintf('Loading session data from: %s\n', data_file);
load(data_file, 'session_data');

%% Step 2: Calculate and Store Diagnostic Metrics (if needed)
if ~isfield(session_data, 'metrics')
giveFeed('Step 2: Calculating and storing diagnostic metrics...');

% Calculate baseline firing rates and add to session_data
session_data.metrics.baseline_frs = calculate_baseline_fr(session_data);

% Calculate waveform metrics and add to session_data
nClusters = height(session_data.spikes.cluster_info.cluster_id);
for i_cluster = 1:nClusters

    % get current unit's multi-channel waveform:
    mean_waveform = session_data.spikes.wfMeans{i_cluster};

    % find channel with max variance:
    [~,max_var_chan] = max(var(mean_waveform,[],2));

    % Note: Assuming a sampling rate of 30000 Hz
    session_data.metrics.wf_metrics(i_cluster, 1) = ...
        calculate_waveform_metrics(mean_waveform(max_var_chan,:), 30000);
end
giveFeed('Diagnostic metrics calculated and stored.');
end

%% Run neuron screening
giveFeed('Step 3: Running neuron screening...');
if strcmp(manifest.screening_status{session_idx}, 'complete')
    giveFeed(['Screening status is ''complete''. ' ...
        'Loading results from session_data...']);

    % Load pre-computed results from the session_data struct
    selected_neurons = session_data.analysis.selected_neurons;
    if isfield(session_data.analysis, 'scSide')
        scSide = session_data.analysis.scSide;
    end

else
    giveFeed('Screening status is ''pending''. Running screening functions...');
    selected_neurons = []; % Initialize empty
    session_data.metadata.unique_id = unique_id;
    if contains(unique_id, 'SNc')
        giveFeed('Session is SNc type. Running screen_da_neurons...');
        selected_neurons = screen_da_neurons(session_data, unique_id);
        giveFeed('DA neuron screening complete.');

    elseif contains(unique_id, 'SC')
        giveFeed('Session is SC type. Running screen_sc_neurons...');

        % screen_sc_neurons will now determine scSide and calculate significance
        [selected_neurons, sig_epoch_comp, scSide] = screen_sc_neurons(session_data);

        % Store results in a new 'analysis' field
        session_data.analysis.scSide = scSide;
        session_data.analysis.sig_epoch_comparison = sig_epoch_comp;

        giveFeed(sprintf('SC neuron screening complete. Determined scSide: %s.', scSide));
    else
        error('Unknown session type. Cannot determine which screening function to run.');
    end
    
    % Add selected_neurons to the analysis field
    session_data.analysis.selected_neurons = selected_neurons;

    % Save the updated session_data struct
    giveFeed('Saving screening results to session_data.mat...');
    save(data_file, 'session_data', '-v7.3');

    % Update manifest status in memory
    manifest.screening_status{session_idx} = 'complete';
    giveFeed('Manifest screening_status updated to ''complete''.');
end


%% Run Core Data Preparation
% This is the main function being tested
giveFeed('Step 4: Preparing core data...');
if strcmp(manifest.dataprep_status{session_idx}, 'complete')
    giveFeed('Data prep status is ''complete''. Loading core_data from session_data...');
    core_data = session_data.analysis.core_data;
else
    giveFeed('Data prep status is ''pending''. Running prepare_core_data...');
    core_data = prepare_core_data(session_data, selected_neurons);

    % Save the core_data struct into the session_data.analysis field
    session_data.analysis.core_data = core_data;

    % Save the updated session_data struct
    giveFeed('Saving core_data to session_data.mat...');
    save(data_file, 'session_data', '-v7.3');

    % Update manifest status in memory
    manifest.dataprep_status{session_idx} = 'complete';
    giveFeed('Manifest dataprep_status updated to ''complete''.');
end

%% Define Task Conditions
% Create logical masks for different trial conditions
giveFeed('Step 5: Defining task conditions...');
[conditions, is_av_session] = define_task_conditions(session_data.trialInfo, ...
    session_data.eventTimes, session_data.metadata.unique_id);

%% Generate Verification Plots
giveFeed('Step 6: Generating verification plots...');
fig = figure('Position', [100, 100, 1200, 700]);
fig.PaperPositionMode = 'auto';

% Define colors for conditions
colors = struct();
colors.normal_dist = [0.2, 0.2, 0.8]; % Blue
colors.uniform_dist = [0.8, 0.2, 0.2]; % Red
colors.norm_common = [0.5, 0.5, 0.5]; % Gray
colors.norm_rare_high = [0.9, 0.6, 0];   % Orange

% Plotting parameters
psth_smoothing_width = 5; % bins

% --- Pre-calculate PSTHs and determine global color scale for heatmaps ---
if any(selected_neurons)
    all_psth_matrices = [];

    % --- Column 1: Aligned to CUE_ON ---
    rates = core_data.spikes.CUE_ON.rates;
    [grand_mean_norm, psth_norm] = calculate_mean_psth(rates, conditions.is_normal_dist);
    [grand_mean_unif, psth_unif] = calculate_mean_psth(rates, conditions.is_uniform_dist);
    all_psth_matrices = [all_psth_matrices; psth_norm; psth_unif];

    % --- Column 2: Aligned to outcomeOn ---
    rates = core_data.spikes.outcomeOn.rates;
    [grand_mean_common, psth_common] = calculate_mean_psth(rates, conditions.is_norm_common);
    [grand_mean_rare, psth_rare] = calculate_mean_psth(rates, conditions.is_norm_rare_high);
    all_psth_matrices = [all_psth_matrices; psth_common; psth_rare];

    % --- Column 3: Aligned to reward ---
    rates = core_data.spikes.reward.rates;
    [grand_mean_common_rew, psth_common_rew] = calculate_mean_psth(rates, conditions.is_norm_common);
    [grand_mean_rare_rew, psth_rare_rew] = calculate_mean_psth(rates, conditions.is_norm_rare_high);
    all_psth_matrices = [all_psth_matrices; psth_common_rew; psth_rare_rew];

    % Determine the common color scale
    max_rate = max(all_psth_matrices, [], 'all');
    c_lims = [0, max_rate];
end

% --- Column 1: Aligned to CUE_ON ---
event_name = 'CUE_ON';
time_vec = core_data.spikes.(event_name).time_vector;
xlim_win = core_data.spikes.(event_name).window;

% Top Row: Per-neuron heatmaps
if any(selected_neurons)
    h(1,1) = mySubPlot([6,3,1]);
    imagesc(time_vec, 1:size(psth_norm, 1), psth_norm);
    clim(c_lims);
    ylabel('Neurons');

    h(2,1) = mySubPlot([6,3,4]);
    imagesc(time_vec, 1:size(psth_unif, 1), psth_unif);
    clim(c_lims);
    ylabel('Neurons');
end

% Middle Row: Grand-average PSTH
h(3,1) = mySubPlot([3,3,4]);
hold on;
if any(selected_neurons)
    barStairsFill(time_vec, grand_mean_norm, 'FaceColor', colors.normal_dist, 'EdgeColor', 'none');
    barStairsFill(time_vec, grand_mean_unif, 'FaceColor', colors.uniform_dist, 'EdgeColor', 'none');
    ylabel('Firing Rate (spikes/s)');
    xlim(xlim_win);
    grid on;
end
xline(0, 'k--');

% Bottom Row: Pupil Trace
h(4,1) = mySubPlot([3,3,7]);
hold on;
traces = core_data.pupil.(event_name).traces;
time_vec_pupil = core_data.pupil.(event_name).time_vector;
xlim_win_pupil = core_data.pupil.(event_name).window;
mean_trace_norm = mean(traces(conditions.is_normal_dist, :), 1, 'omitnan');
plot(time_vec_pupil, mean_trace_norm, 'Color', colors.normal_dist, 'LineWidth', 1.5);
mean_trace_unif = mean(traces(conditions.is_uniform_dist, :), 1, 'omitnan');
plot(time_vec_pupil, mean_trace_unif, 'Color', colors.uniform_dist, 'LineWidth', 1.5);
ylabel('Pupil Size (norm.)');
xlabel(['Time from ' strrep(event_name, '_', ' ') ' (s)']);
xlim(xlim_win_pupil);
grid on;
xline(0, 'k--');

% --- Column 2: Aligned to outcomeOn ---
event_name = 'outcomeOn';
time_vec = core_data.spikes.(event_name).time_vector;
xlim_win = core_data.spikes.(event_name).window;

% Top Row: Per-neuron heatmaps
if any(selected_neurons)
    h(1,2) = mySubPlot([6,3,2]);
    imagesc(time_vec, 1:size(psth_common, 1), psth_common);
    clim(c_lims);

    h(2,2) = mySubPlot([6,3,5]);
    imagesc(time_vec, 1:size(psth_rare, 1), psth_rare);
    clim(c_lims);
end

% Middle Row: Grand-average PSTH
h(3,2) = mySubPlot([3,3,5]);
hold on;
if any(selected_neurons)
    barStairsFill(time_vec, grand_mean_common, 'FaceColor', colors.norm_common, 'EdgeColor', 'none');
    barStairsFill(time_vec, grand_mean_rare, 'FaceColor', colors.norm_rare_high, 'EdgeColor', 'none');
    xlim(xlim_win);
    grid on;
end
xline(0, 'k--');

% Bottom Row: Pupil Trace
h(4,2) = mySubPlot([3,3,8]);
hold on;
traces = core_data.pupil.(event_name).traces;
time_vec_pupil = core_data.pupil.(event_name).time_vector;
xlim_win_pupil = core_data.pupil.(event_name).window;
mean_trace_common = mean(traces(conditions.is_norm_common, :), 1, 'omitnan');
plot(time_vec_pupil, mean_trace_common, 'Color', colors.norm_common, 'LineWidth', 1.5);
mean_trace_rare = mean(traces(conditions.is_norm_rare_high, :), 1, 'omitnan');
plot(time_vec_pupil, mean_trace_rare, 'Color', colors.norm_rare_high, 'LineWidth', 1.5);
xlabel(['Time from ' strrep(event_name, '_', ' ') ' (s)']);
xlim(xlim_win_pupil);
grid on;
xline(0, 'k--');

% --- Column 3: Aligned to reward ---
event_name = 'reward';
time_vec = core_data.spikes.(event_name).time_vector;
xlim_win = core_data.spikes.(event_name).window;

% Top Row: Per-neuron heatmaps
if any(selected_neurons)
    h(1,3) = mySubPlot([6,3,3]);
    imagesc(time_vec, 1:size(psth_common_rew, 1), psth_common_rew);
    clim(c_lims);

    h(2,3) = mySubPlot([6,3,6]);
    imagesc(time_vec, 1:size(psth_rare_rew, 1), psth_rare_rew);
    clim(c_lims);
    colorbar;
end

% Middle Row: Grand-average PSTH
h(3,3) = mySubPlot([3,3,6]);
hold on;
if any(selected_neurons)
    barStairsFill(time_vec, grand_mean_common_rew, 'FaceColor', colors.norm_common, 'EdgeColor', 'none');
    barStairsFill(time_vec, grand_mean_rare_rew, 'FaceColor', colors.norm_rare_high, 'EdgeColor', 'none');
    xlim(xlim_win);
    grid on;
end
xline(0, 'k--');

% Bottom Row: Pupil Trace
h(4,3) = mySubPlot([3,3,9]);
hold on;
traces = core_data.pupil.(event_name).traces;
time_vec_pupil = core_data.pupil.(event_name).time_vector;
xlim_win_pupil = core_data.pupil.(event_name).window;
mean_trace_common = mean(traces(conditions.is_norm_common, :), 1, 'omitnan');
plot(time_vec_pupil, mean_trace_common, 'Color', colors.norm_common, 'LineWidth', 1.5);
mean_trace_rare = mean(traces(conditions.is_norm_rare_high, :), 1, 'omitnan');
plot(time_vec_pupil, mean_trace_rare, 'Color', colors.norm_rare_high, 'LineWidth', 1.5);
xlabel(['Time from ' strrep(event_name, '_', ' ') ' (s)']);
xlim(xlim_win_pupil);
grid on;
xline(0, 'k--');

% --- De-clutter axes for clarity ---
if any(selected_neurons)
    % Remove X-tick labels from all but the bottom row (pupil)
    all_axes = findobj(fig, 'Type', 'axes');
    psth_and_heatmap_axes = all_axes(~ismember(all_axes, h(4,:)));
    set(psth_and_heatmap_axes, 'XTickLabel', []);

    % Remove Y-tick labels from all but the left-most column
    set(h(:, 2:3), 'YTickLabel', []);
end

sgtitle(sprintf('Core Data Verification: %s', unique_id), 'Interpreter', 'none');
giveFeed('Done.');

%% Generate SPE-Focused Diagnostic Figure (if applicable)
if is_av_session
    giveFeed('Step 6b: Generating SPE-focused diagnostic plots...');
    fig2 = figure('Position', [100, 100, 1200, 700]);
    fig2.PaperPositionMode = 'auto';

    % Define colors for SPE conditions
    colors.flicker_surprising = [0.8, 0.2, 0.8]; % Magenta
    colors.flicker_certain = [0.2, 0.8, 0.2]; % Green

    % --- Column 1: Aligned to CUE_ON ---
    event_name = 'CUE_ON';
    time_vec = core_data.spikes.(event_name).time_vector;
    xlim_win = core_data.spikes.(event_name).window;

    % Top Row: Neuronal PSTH
    h2(1) = axes('Position', pos{1,1});
    hold on;
    if any(selected_neurons)
        rates = core_data.spikes.(event_name).rates;
        mean_psth_surprising = mean(squeeze(mean(rates(:, ...
            conditions.is_flicker_surprising, :), 1, 'omitnan')), 1, ...
            'omitnan');
        hOut = barStairsFill(time_vec, zeros(size(mean_psth_surprising)), ...
            mean_psth_surprising);
        delete(hOut(1:2));
        set(hOut(3), 'Color', colors.flicker_surprising)
        mean_psth_certain = mean(squeeze(mean(rates(:, ...
            conditions.is_flicker_certain, :), 1, 'omitnan')), 1, ...
            'omitnan');
        hOut = barStairsFill(time_vec, ...
            zeros(size(mean_psth_certain)), mean_psth_certain);
        delete(hOut(1:2));
        set(hOut(3), 'Color', colors.flicker_certain)
        ylabel('Firing Rate (spikes/s)');
        xlim(xlim_win);
        grid on;
        legend({'Flicker Surprising', 'Flicker Certain'}, ...
            'Location', 'northeast');
    end
    xline(0, 'k--');

    % Bottom Row: Pupil Trace
    h2(4) = axes('Position', pos{2,1});
    hold on;
    traces = core_data.pupil.(event_name).traces;
    time_vec = core_data.pupil.(event_name).time_vector;
    xlim_win = core_data.pupil.(event_name).window;
    mean_trace_surprising = mean(traces(conditions.is_flicker_surprising, :), 1, 'omitnan');
    plot(time_vec, mean_trace_surprising, 'Color', colors.flicker_surprising, 'LineWidth', 1.5);
    mean_trace_certain = mean(traces(conditions.is_flicker_certain, :), 1, 'omitnan');
    plot(time_vec, mean_trace_certain, 'Color', colors.flicker_certain, 'LineWidth', 1.5);
    ylabel('Normalized Pupil Size');
    xlabel('Time from Cue On (s)');
    xlim(xlim_win);
    grid on;
    legend({'Flicker Surprising', 'Flicker Certain'}, 'Location', 'northeast');
    xline(0, 'k--');

    % --- Column 2: Aligned to outcomeOn ---
    event_name = 'outcomeOn';
    time_vec = core_data.spikes.(event_name).time_vector;
    xlim_win = core_data.spikes.(event_name).window;

    % Top Row: Neuronal PSTH
    h2(2) = axes('Position', pos{1,2});
    hold on;
    if any(selected_neurons)
        rates = core_data.spikes.(event_name).rates;
        mean_psth_surprising = mean(squeeze(mean(rates(:, ...
            conditions.is_flicker_surprising, :), 1, 'omitnan')), 1, ...
            'omitnan');
        hOut = barStairsFill(time_vec, zeros(size(mean_psth_surprising)), ...
            mean_psth_surprising);
        delete(hOut(1:2));
        set(hOut(3), 'Color', colors.flicker_surprising)
        mean_psth_certain = mean(squeeze(mean(rates(:, ...
            conditions.is_flicker_certain, :), 1, 'omitnan')), 1, ...
            'omitnan');
        hOut = barStairsFill(time_vec, zeros(size(mean_psth_certain)), ...
            mean_psth_certain);
        delete(hOut(1:2));
        set(hOut(3), 'Color', colors.flicker_certain)
        xlim(xlim_win);
        grid on;
    end
    xline(0, 'k--');

    % Bottom Row: Pupil Trace
    h2(5) = axes('Position', pos{2,2});
    hold on;
    traces = core_data.pupil.(event_name).traces;
    time_vec = core_data.pupil.(event_name).time_vector;
    xlim_win = core_data.pupil.(event_name).window;
    mean_trace_surprising = mean(traces(conditions.is_flicker_surprising, :), 1, 'omitnan');
    plot(time_vec, mean_trace_surprising, 'Color', colors.flicker_surprising, 'LineWidth', 1.5);
    mean_trace_certain = mean(traces(conditions.is_flicker_certain, :), 1, 'omitnan');
    plot(time_vec, mean_trace_certain, 'Color', colors.flicker_certain, 'LineWidth', 1.5);
    xlabel('Time from Outcome On (s)');
    xlim(xlim_win);
    grid on;
    xline(0, 'k--');

    % --- Column 3: Aligned to reward ---
    event_name = 'reward';
    time_vec = core_data.spikes.(event_name).time_vector;
    xlim_win = core_data.spikes.(event_name).window;

    % Top Row: Neuronal PSTH
    h2(3) = axes('Position', pos{1,3});
    hold on;
    if any(selected_neurons)
        rates = core_data.spikes.(event_name).rates;
        mean_psth_surprising = mean(squeeze(mean(rates(:, ...
            conditions.is_flicker_surprising, :), 1, 'omitnan')), 1, ...
            'omitnan');
        hOut = barStairsFill(time_vec, zeros(size(mean_psth_surprising)), ...
            mean_psth_surprising);
        delete(hOut(1:2));
        set(hOut(3), 'Color', colors.flicker_surprising)
        mean_psth_certain = mean(squeeze(mean(rates(:, ...
            conditions.is_flicker_certain, :), 1, 'omitnan')), 1, ...
            'omitnan');
        hOut = barStairsFill(time_vec, zeros(size(mean_psth_certain)), ...
            mean_psth_certain);
        delete(hOut(1:2));
        set(hOut(3), 'Color', colors.flicker_certain)
        xlim(xlim_win);
        grid on;
        legend({'Flicker Surprising', 'Flicker Certain'}, ...
            'Location', 'northeast');
    end
    xline(0, 'k--');

    % Bottom Row: Pupil Trace
    h2(6) = axes('Position', pos{2,3});
    hold on;
    traces = core_data.pupil.(event_name).traces;
    time_vec = core_data.pupil.(event_name).time_vector;
    xlim_win = core_data.pupil.(event_name).window;
    mean_trace_surprising = mean(traces(conditions.is_flicker_surprising, :), 1, 'omitnan');
    plot(time_vec, mean_trace_surprising, 'Color', colors.flicker_surprising, 'LineWidth', 1.5);
    mean_trace_certain = mean(traces(conditions.is_flicker_certain, :), 1, 'omitnan');
    plot(time_vec, mean_trace_certain, 'Color', colors.flicker_certain, 'LineWidth', 1.5);
    xlabel('Time from Reward (s)');
    xlim(xlim_win);
    grid on;
    legend({'Flicker Surprising', 'Flicker Certain'}, 'Location', 'northeast');
    xline(0, 'k--');

    % --- De-clutter axes ---
    set(h2(1:3), 'XTickLabel', []);
    set(h2([2,3,5,6]), 'YTickLabel', []);

    sgtitle(sprintf('SPE-Focused Data Verification: %s', unique_id), 'Interpreter', 'none');
    giveFeed('Done with SPE-focused plots.');

    % --- Link Y-axes across both figures for direct comparison ---
    giveFeed('Step 6c: Linking Y-axes across figures...');

    % Link PSTH plots (h from fig 1, h2 from fig 2)
    psth_axes = [h(1:3), h2(1:3)];
    [~, y_lims_psth] = outerLims(psth_axes);
    set(psth_axes, 'YLim', y_lims_psth);

    % Link Pupil plots (h from fig 1, h2 from fig 2)
    pupil_axes = [h(4:6), h2(4:6)];
    [~, y_lims_pupil] = outerLims(pupil_axes);
    set(pupil_axes, 'YLim', y_lims_pupil);
end

%% Finalize and Save Manifest
giveFeed('Step 7: Saving updated manifest...');
writetable(manifest, manifest_path);
giveFeed('Manifest saved. Script finished.');

%% Helper function
function [grand_mean_psth, per_neuron_psth] = calculate_mean_psth(rates, condition_mask)
%calculate_mean_psth Calculates the mean PSTH, excluding silent trials.
%   This function calculates the grand-average PSTH across all neurons and
%   the average PSTH for each individual neuron. It filters out trials
%   where a neuron had a total spike count of zero across the analysis
%   window on a per-neuron basis.
%
%   Args:
%       rates: An [nNeurons x nTrials x nTimeBins] matrix of spike rates.
%       condition_mask: A logical vector indicating which trials to include.
%
%   Returns:
%       grand_mean_psth: A [1 x nTimeBins] vector of the grand-average PSTH.
%       per_neuron_psth: An [nNeurons x nTimeBins] matrix of per-neuron PSTHs.

    n_neurons = size(rates, 1);
    n_time_bins = size(rates, 3);

    per_neuron_psth = zeros(n_neurons, n_time_bins);

    for i_neuron = 1:n_neurons
        % Get the rates for the current neuron and the specified condition
        neuron_rates = squeeze(rates(i_neuron, condition_mask, :));

        % Find trials where the neuron was not silent
        active_trials = sum(neuron_rates, 2) > 0;

        % Calculate the mean PSTH for the active trials
        if any(active_trials)
            per_neuron_psth(i_neuron, :) = mean(neuron_rates(active_trials, :), 1);
        else
            % If the neuron is silent in all trials, its mean PSTH is zero
            per_neuron_psth(i_neuron, :) = zeros(1, n_time_bins);
        end
    end

    % Calculate the grand-average PSTH from the per-neuron averages
    grand_mean_psth = mean(per_neuron_psth, 1);
end
