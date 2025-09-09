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
alignment_events = fieldnames(core_data.spikes); % {'CUE_ON', 'outcomeOn', 'reward'}

%% Generate Verification Plots
fprintf('Generating verification plots...\n');
fig = figure('Position', [100, 100, 1500, 800]);
fig.PaperPositionMode = 'auto';

% Plotting parameters
psth_smoothing_width = 5; % bins

for i_event = 1:numel(alignment_events)
    event_name = alignment_events{i_event};

    % --- Top Row: Neuronal PSTHs ---
    mySubPlot(2, 3, i_event, 'shape', 'wide');

    if any(selected_neurons)
        % Extract data for the current event
        rates = core_data.spikes.(event_name).rates; % (neurons x trials x time)
        time_vec = core_data.spikes.(event_name).time_vector;

        % Calculate the grand average PSTH across all neurons and trials
        % Note: Taking the mean twice (first over neurons, then over trials)
        % is equivalent to a grand average but handles NaNs better.
        mean_psth = mean(squeeze(mean(rates, 1, 'omitnan')), 1, 'omitnan');

        % Smooth the PSTH for visualization
        smoothed_psth = movmean(mean_psth, psth_smoothing_width);

        % Plotting
        plot(time_vec, smoothed_psth, 'k', 'LineWidth', 2);
        hold on;
        xline(0, 'r--', 'LineWidth', 1); % Mark alignment time

        % Formatting
        title(['PSTH: ' strrep(event_name, '_', ' ')]);
        xlabel('Time from event (s)');
        ylabel('Firing Rate (spikes/s)');
        xlim(core_data.spikes.(event_name).window);
        grid on;
    else
        title(['PSTH: ' strrep(event_name, '_', ' ') ' (No Neurons)']);
        axis off;
    end

    % --- Bottom Row: Pupil Traces ---
    mySubPlot(2, 3, i_event + 3, 'shape', 'wide');

    % Extract data for the current event
    traces = core_data.pupil.(event_name).traces; % (trials x time)
    time_vec = core_data.pupil.(event_name).time_vector;

    % Calculate the grand average pupil trace across all trials
    mean_trace = mean(traces, 1, 'omitnan');

    % Plotting
    plot(time_vec, mean_trace, 'b', 'LineWidth', 2);
    hold on;
    xline(0, 'r--', 'LineWidth', 1); % Mark alignment time

    % Formatting
    title(['Pupil Trace: ' strrep(event_name, '_', ' ')]);
    xlabel('Time from event (s)');
    ylabel('Normalized Pupil Size');
    xlim(core_data.pupil.(event_name).window);
    grid on;
end

sgtitle('Core Data Preparation Verification');
fprintf('Done.\n');

end
