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

% Add utility functions to the MATLAB path
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Load Data
% Define the project root and construct the path to a sample data file.
% NOTE: The user should replace 'sample_session.mat' with a real session
% file from their 'data/raw' directory.
project_root = fullfile(findOneDrive, 'Code', 'tokens-analysis-pipeline');
data_file = fullfile(project_root, 'data', 'raw', 'sample_session.mat');

% Check if the sample file exists
if ~exist(data_file, 'file')
    error('test_core_data_preparation:sampleNotFound', ...
          ['Sample data file not found. Please specify a valid path to a \n' ...
           'session_data.mat file in the "Load Data" section of this script.']);
end

fprintf('Loading session data from: %s\n', data_file);
load(data_file, 'session_data');

%% Run Neuron Screening
% Use the standard screening function to select neurons for analysis
fprintf('Screening SC neurons...\n');
[selected_neurons, ~, ~] = screen_sc_neurons(session_data);

% Check if any neurons were selected
if ~any(selected_neurons)
    warning('No neurons passed the screening criteria. Plots will be empty.');
end

%% Run Core Data Preparation
% This is the main function being tested
fprintf('Running core data preparation...\n');
core_data = prepare_core_data(session_data, selected_neurons);

%% Define Task Conditions
% Create logical masks for different trial conditions
fprintf('Defining task conditions...\n');
conditions = define_task_conditions(session_data.trialInfo, ...
    session_data.eventTimes);

%% Generate Verification Plots
fprintf('Generating verification plots...\n');
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

% --- Create all subplots first ---
h(1) = mySubPlot([2, 3, 1]);
h(4) = mySubPlot([2, 3, 4]);
h(2) = mySubPlot([2, 3, 2]);
h(5) = mySubPlot([2, 3, 5]);
h(3) = mySubPlot([2, 3, 3]);
h(6) = mySubPlot([2, 3, 6]);

% --- Column 1: Aligned to CUE_ON ---
event_name = 'CUE_ON';
time_vec = core_data.spikes.(event_name).time_vector;
xlim_win = core_data.spikes.(event_name).window;

% Top Row: Neuronal PSTH
axes(h(1));
hold on;
if any(selected_neurons)
    rates = core_data.spikes.(event_name).rates;

    % Normal distribution
    mean_psth_norm = mean(squeeze(mean(rates(:, conditions.is_normal_dist, :), 1, 'omitnan')), 1, 'omitnan');
    plot(time_vec, movmean(mean_psth_norm, psth_smoothing_width), 'Color', colors.normal_dist, 'LineWidth', 1.5);

    % Uniform distribution
    mean_psth_unif = mean(squeeze(mean(rates(:, conditions.is_uniform_dist, :), 1, 'omitnan')), 1, 'omitnan');
    plot(time_vec, movmean(mean_psth_unif, psth_smoothing_width), 'Color', colors.uniform_dist, 'LineWidth', 1.5);

    ylabel('Firing Rate (spikes/s)');
    xlim(xlim_win);
    grid on;
    legend({'Normal Dist', 'Uniform Dist'}, 'Location', 'northeast');
end
title('PSTH: Cue On');
xline(0, 'k--');

% Bottom Row: Pupil Trace
axes(h(4));
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
title('Pupil: Cue On');
xline(0, 'k--');

% --- Column 2: Aligned to outcomeOn ---
event_name = 'outcomeOn';
time_vec = core_data.spikes.(event_name).time_vector;
xlim_win = core_data.spikes.(event_name).window;

% Top Row: Neuronal PSTH
axes(h(2));
hold on;
if any(selected_neurons)
    rates = core_data.spikes.(event_name).rates;

    % Common reward
    mean_psth_common = mean(squeeze(mean(rates(:, conditions.is_norm_common, :), 1, 'omitnan')), 1, 'omitnan');
    plot(time_vec, movmean(mean_psth_common, psth_smoothing_width), 'Color', colors.norm_common, 'LineWidth', 1.5);

    % Rare high reward
    mean_psth_rare = mean(squeeze(mean(rates(:, conditions.is_norm_rare_high, :), 1, 'omitnan')), 1, 'omitnan');
    plot(time_vec, movmean(mean_psth_rare, psth_smoothing_width), 'Color', colors.norm_rare_high, 'LineWidth', 1.5);

    xlim(xlim_win);
    grid on;
end
title('PSTH: Outcome On');
xline(0, 'k--');

% Bottom Row: Pupil Trace
axes(h(5));
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
title('Pupil: Outcome On');
xline(0, 'k--');

% --- Column 3: Aligned to reward ---
event_name = 'reward';
time_vec = core_data.spikes.(event_name).time_vector;
xlim_win = core_data.spikes.(event_name).window;

% Top Row: Neuronal PSTH
axes(h(3));
hold on;
if any(selected_neurons)
    rates = core_data.spikes.(event_name).rates;

    % Common reward
    mean_psth_common = mean(squeeze(mean(rates(:, conditions.is_norm_common, :), 1, 'omitnan')), 1, 'omitnan');
    plot(time_vec, movmean(mean_psth_common, psth_smoothing_width), 'Color', colors.norm_common, 'LineWidth', 1.5);

    % Rare high reward
    mean_psth_rare = mean(squeeze(mean(rates(:, conditions.is_norm_rare_high, :), 1, 'omitnan')), 1, 'omitnan');
    plot(time_vec, movmean(mean_psth_rare, psth_smoothing_width), 'Color', colors.norm_rare_high, 'LineWidth', 1.5);

    xlim(xlim_win);
    grid on;
    legend({'Common Reward', 'Rare High Reward'}, 'Location', 'northeast');
end
title('PSTH: Reward');
xline(0, 'k--');

% Bottom Row: Pupil Trace
axes(h(6));
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
title('Pupil: Reward');
xline(0, 'k--');
legend({'Common Reward', 'Rare High Reward'}, 'Location', 'northeast');

sgtitle('Core Data Preparation Verification: Condition Comparisons');
fprintf('Done.\n');
