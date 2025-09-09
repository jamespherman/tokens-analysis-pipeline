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
unique_id = 'Feynman_08_05_2025_SNc'; % Hardcoded session
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

% --- Manually define subplot positions for proportional widths ---
time_windows = [2, 2, 5.5]; % Durations for Cue, Outcome, Reward
total_time = sum(time_windows);
proportions = time_windows / total_time;

% Layout settings
left_margin = 0.07;
right_margin = 0.03;
top_margin = 0.1;
bottom_margin = 0.1;
h_gap = 0.04; % Horizontal gap
v_gap = 0.08; % Vertical gap

plot_area_width = 1 - left_margin - right_margin - 2 * h_gap;
plot_area_height = 1 - top_margin - bottom_margin - v_gap;

col_widths = plot_area_width * proportions;
row_height = plot_area_height / 2;

% Calculate positions [left, bottom, width, height]
pos{1,1} = [left_margin, bottom_margin + row_height + v_gap, col_widths(1), row_height];
pos{2,1} = [left_margin, bottom_margin, col_widths(1), row_height];

pos{1,2} = [left_margin + col_widths(1) + h_gap, bottom_margin + row_height + v_gap, col_widths(2), row_height];
pos{2,2} = [left_margin + col_widths(1) + h_gap, bottom_margin, col_widths(2), row_height];

pos{1,3} = [left_margin + col_widths(1) + col_widths(2) + 2*h_gap, bottom_margin + row_height + v_gap, col_widths(3), row_height];
pos{2,3} = [left_margin + col_widths(1) + col_widths(2) + 2*h_gap, bottom_margin, col_widths(3), row_height];


% --- Column 1: Aligned to CUE_ON ---
event_name = 'CUE_ON';
time_vec = core_data.spikes.(event_name).time_vector;
xlim_win = core_data.spikes.(event_name).window;

% Top Row: Neuronal PSTH
h(1) = axes('Position', pos{1,1});
hold on;
if any(selected_neurons)
    rates = core_data.spikes.(event_name).rates;

    % Normal distribution
    mean_psth_norm = mean(squeeze(mean(rates(:, conditions.is_normal_dist, :), 1, 'omitnan')), 1, 'omitnan');
    barStairsFill(time_vec, movmean(mean_psth_norm, psth_smoothing_width), 'FaceColor', colors.normal_dist, 'EdgeColor', 'none', 'FaceAlpha', 0.7);

    % Uniform distribution
    mean_psth_unif = mean(squeeze(mean(rates(:, conditions.is_uniform_dist, :), 1, 'omitnan')), 1, 'omitnan');
    barStairsFill(time_vec, movmean(mean_psth_unif, psth_smoothing_width), 'FaceColor', colors.uniform_dist, 'EdgeColor', 'none', 'FaceAlpha', 0.7);

    ylabel('Firing Rate (spikes/s)');
    xlim(xlim_win);
    grid on;
    legend({'Normal Dist', 'Uniform Dist'}, 'Location', 'northeast');
end
xline(0, 'k--');

% Bottom Row: Pupil Trace
h(4) = axes('Position', pos{2,1});
hold on;
traces = core_data.pupil.(event_name).traces;
time_vec = core_data.pupil.(event_name).time_vector;
xlim_win = core_data.pupil.(event_name).window;

% Normal distribution
mean_trace_norm = mean(traces(conditions.is_normal_dist, :), 1, 'omitnan');
plot(time_vec, mean_trace_norm, 'Color', colors.normal_dist, 'LineWidth', 1.5);

% Uniform distribution
mean_trace_unif = mean(traces(conditions.is_uniform_dist, :), 1, 'omitnan');
plot(time_vec, mean_trace_unif, 'Color', colors.uniform_dist, 'LineWidth', 1.5);

ylabel('Normalized Pupil Size');
xlabel('Time from Cue On (s)');
xlim(xlim_win);
grid on;
legend({'Normal Dist', 'Uniform Dist'}, 'Location', 'northeast');
xline(0, 'k--');

% --- Column 2: Aligned to outcomeOn ---
event_name = 'outcomeOn';
time_vec = core_data.spikes.(event_name).time_vector;
xlim_win = core_data.spikes.(event_name).window;

% Top Row: Neuronal PSTH
h(2) = axes('Position', pos{1,2});
hold on;
if any(selected_neurons)
    rates = core_data.spikes.(event_name).rates;

    % Common reward
    mean_psth_common = mean(squeeze(mean(rates(:, conditions.is_norm_common, :), 1, 'omitnan')), 1, 'omitnan');
    barStairsFill(time_vec, movmean(mean_psth_common, psth_smoothing_width), 'FaceColor', colors.norm_common, 'EdgeColor', 'none', 'FaceAlpha', 0.7);

    % Rare high reward
    mean_psth_rare = mean(squeeze(mean(rates(:, conditions.is_norm_rare_high, :), 1, 'omitnan')), 1, 'omitnan');
    barStairsFill(time_vec, movmean(mean_psth_rare, psth_smoothing_width), 'FaceColor', colors.norm_rare_high, 'EdgeColor', 'none', 'FaceAlpha', 0.7);

    xlim(xlim_win);
    grid on;
end
xline(0, 'k--');

% Bottom Row: Pupil Trace
h(5) = axes('Position', pos{2,2});
hold on;
traces = core_data.pupil.(event_name).traces;
time_vec = core_data.pupil.(event_name).time_vector;
xlim_win = core_data.pupil.(event_name).window;

% Common reward
mean_trace_common = mean(traces(conditions.is_norm_common, :), 1, 'omitnan');
plot(time_vec, mean_trace_common, 'Color', colors.norm_common, 'LineWidth', 1.5);

% Rare high reward
mean_trace_rare = mean(traces(conditions.is_norm_rare_high, :), 1, 'omitnan');
plot(time_vec, mean_trace_rare, 'Color', colors.norm_rare_high, 'LineWidth', 1.5);

xlabel('Time from Outcome On (s)');
xlim(xlim_win);
grid on;
xline(0, 'k--');

% --- Column 3: Aligned to reward ---
event_name = 'reward';
time_vec = core_data.spikes.(event_name).time_vector;
xlim_win = core_data.spikes.(event_name).window;

% Top Row: Neuronal PSTH
h(3) = axes('Position', pos{1,3});
hold on;
if any(selected_neurons)
    rates = core_data.spikes.(event_name).rates;

    % Common reward
    mean_psth_common = mean(squeeze(mean(rates(:, conditions.is_norm_common, :), 1, 'omitnan')), 1, 'omitnan');
    barStairsFill(time_vec, movmean(mean_psth_common, psth_smoothing_width), 'FaceColor', colors.norm_common, 'EdgeColor', 'none', 'FaceAlpha', 0.7);

    % Rare high reward
    mean_psth_rare = mean(squeeze(mean(rates(:, conditions.is_norm_rare_high, :), 1, 'omitnan')), 1, 'omitnan');
    barStairsFill(time_vec, movmean(mean_psth_rare, psth_smoothing_width), 'FaceColor', colors.norm_rare_high, 'EdgeColor', 'none', 'FaceAlpha', 0.7);

    xlim(xlim_win);
    grid on;
    legend({'Common Reward', 'Rare High Reward'}, 'Location', 'northeast');
end
xline(0, 'k--');

% Bottom Row: Pupil Trace
h(6) = axes('Position', pos{2,3});
hold on;
traces = core_data.pupil.(event_name).traces;
time_vec = core_data.pupil.(event_name).time_vector;
xlim_win = core_data.pupil.(event_name).window;

% Common reward
mean_trace_common = mean(traces(conditions.is_norm_common, :), 1, 'omitnan');
plot(time_vec, mean_trace_common, 'Color', colors.norm_common, 'LineWidth', 1.5);

% Rare high reward
mean_trace_rare = mean(traces(conditions.is_norm_rare_high, :), 1, 'omitnan');
plot(time_vec, mean_trace_rare, 'Color', colors.norm_rare_high, 'LineWidth', 1.5);

xlabel('Time from Reward (s)');
xlim(xlim_win);
grid on;
legend({'Common Reward', 'Rare High Reward'}, 'Location', 'northeast');
xline(0, 'k--');

% --- De-clutter axes ---
set(h(1:3), 'XTickLabel', []);
set(h([2,3,5,6]), 'YTickLabel', []);

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
        mean_psth_surprising = mean(squeeze(mean(rates(:, conditions.is_flicker_surprising, :), 1, 'omitnan')), 1, 'omitnan');
        barStairsFill(time_vec, movmean(mean_psth_surprising, psth_smoothing_width), 'FaceColor', colors.flicker_surprising, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
        mean_psth_certain = mean(squeeze(mean(rates(:, conditions.is_flicker_certain, :), 1, 'omitnan')), 1, 'omitnan');
        barStairsFill(time_vec, movmean(mean_psth_certain, psth_smoothing_width), 'FaceColor', colors.flicker_certain, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
        ylabel('Firing Rate (spikes/s)');
        xlim(xlim_win);
        grid on;
        legend({'Flicker Surprising', 'Flicker Certain'}, 'Location', 'northeast');
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
        mean_psth_surprising = mean(squeeze(mean(rates(:, conditions.is_flicker_surprising, :), 1, 'omitnan')), 1, 'omitnan');
        barStairsFill(time_vec, movmean(mean_psth_surprising, psth_smoothing_width), 'FaceColor', colors.flicker_surprising, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
        mean_psth_certain = mean(squeeze(mean(rates(:, conditions.is_flicker_certain, :), 1, 'omitnan')), 1, 'omitnan');
        barStairsFill(time_vec, movmean(mean_psth_certain, psth_smoothing_width), 'FaceColor', colors.flicker_certain, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
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
        mean_psth_surprising = mean(squeeze(mean(rates(:, conditions.is_flicker_surprising, :), 1, 'omitnan')), 1, 'omitnan');
        barStairsFill(time_vec, movmean(mean_psth_surprising, psth_smoothing_width), 'FaceColor', colors.flicker_surprising, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
        mean_psth_certain = mean(squeeze(mean(rates(:, conditions.is_flicker_certain, :), 1, 'omitnan')), 1, 'omitnan');
        barStairsFill(time_vec, movmean(mean_psth_certain, psth_smoothing_width), 'FaceColor', colors.flicker_certain, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
        xlim(xlim_win);
        grid on;
        legend({'Flicker Surprising', 'Flicker Certain'}, 'Location', 'northeast');
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
    y_lims_psth = outerLims(psth_axes, 'y');
    set(psth_axes, 'YLim', y_lims_psth);

    % Link Pupil plots (h from fig 1, h2 from fig 2)
    pupil_axes = [h(4:6), h2(4:6)];
    y_lims_pupil = outerLims(pupil_axes, 'y');
    set(pupil_axes, 'YLim', y_lims_pupil);
end

%% Finalize and Save Manifest
giveFeed('Step 7: Saving updated manifest...');
writetable(manifest, manifest_path);
giveFeed('Manifest saved. Script finished.');
