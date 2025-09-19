%% analyze_baseline_comparison.m
%
%   Compares baseline firing rates to post-event firing rates for a specified
%   condition and alignment events using ROC analysis.
%
% INPUTS:
%   core_data  - A struct containing processed neuronal data.
%   conditions - A struct with logical masks for different trial conditions.
%   is_av_session - A boolean indicating if the session is an AV session.
%
% NAME-VALUE PAIR INPUTS:
%   'condition' - A string specifying the condition to analyze. This is a
%                 required parameter.
%
% OUTPUT:
%   analysis_results - A struct containing the significance matrix for the
%                      specified condition and each alignment.
%
% Author: Jules
% Date: 2025-09-19

function analysis_results = analyze_baseline_comparison(core_data, conditions, is_av_session, varargin)

%% Setup Paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Parse Optional Arguments
p = inputParser;
% The 'condition' parameter is now required. The validation function ensures
% it is a non-empty character string.
addParameter(p, 'condition', '', @(x) ischar(x) && ~isempty(x));
parse(p, varargin{:});

% Error if the 'condition' parameter is not provided by the user.
if any(strcmp(p.UsingDefaults, 'condition'))
    error('The required name-value pair ''condition'' was not provided.');
end
cond_name = p.Results.condition;

%% Analysis Parameters
baseline_window = [-0.2, 0]; % -200ms to 0ms relative to event
roc_perms = 200;

% Define the alignment events
align_events = {'CUE_ON', 'outcomeOn', 'reward'};

n_neurons = size(core_data.spikes.CUE_ON.rates, 1);

%% Pre-calculate Total Workload
total_workload = 0;
total_clusters_in_session = 0;
for i_align_dry = 1:numel(align_events)
    align_name_dry = align_events{i_align_dry};
    time_vector_dry = core_data.spikes.(align_name_dry).time_vector;

    baseline_bins_dry = time_vector_dry >= baseline_window(1) & ...
        time_vector_dry < baseline_window(2);
    post_event_bins_dry = time_vector_dry >= 0;

    n_cols = sum(post_event_bins_dry);

    trial_indices_dry = conditions.(cond_name);
    n_rows = sum(trial_indices_dry);

    % The analysis will skip if there are fewer than 5 trials, so we
    % should not include it in the workload calculation.
    if n_rows >= 5 && (n_rows * sum(baseline_bins_dry)) > 0
        workload_per_neuron = n_rows * n_cols;
        total_workload = total_workload + (workload_per_neuron * ...
            n_neurons);
        total_clusters_in_session = total_clusters_in_session + ...
        n_neurons;
    end
end

%% Initialize Tracking Variables
workload_completed = 0;
clusters_processed = 0;
master_timer = tic;
if total_clusters_in_session > 0
    fprintf('Starting analysis for %d clusters...\n', ...
        total_clusters_in_session);
end

%% Main Analysis Loop
% Iterate through each alignment event
for i_align = 1:numel(align_events)
    align_name = align_events{i_align};
    time_vector = core_data.spikes.(align_name).time_vector;

    % Find time bins corresponding to the baseline and post-event periods
    baseline_bins = time_vector >= baseline_window(1) & ...
        time_vector < baseline_window(2);
    post_event_bins = time_vector >= 0;

    trial_indices = conditions.(cond_name);

    % Get the full rates matrix for this alignment
    rates_all_neurons = core_data.spikes.(align_name).rates;

    % Preallocate storage for significance results for this condition
    n_post_bins = sum(post_event_bins);
    sig_results = nan(n_neurons, n_post_bins);

    % Iterate through each neuron
    for i_neuron = 1:n_neurons
        % Extract rates for the current neuron and condition
        % Resulting matrix is [1 x nTrials x nTimeBins]
        neuron_rates_all_bins = rates_all_neurons(i_neuron, ...
            trial_indices, :);

        % Squeeze to [nTrials x nTimeBins]
        neuron_rates_all_bins = squeeze(neuron_rates_all_bins);

        % Filter out silent trials (all zeros or NaNs)
        is_silent_trial = all(neuron_rates_all_bins == 0 | ...
            isnan(neuron_rates_all_bins), 2);
        active_trials = ~is_silent_trial;
        neuron_rates_all_bins = neuron_rates_all_bins(...
            active_trials, :);

        % If no active trials, skip to the next neuron
        if isempty(neuron_rates_all_bins)
            continue;
        end

        % Construct the baseline distribution (x)
        % Pool values from all baseline bins across all trials
        baseline_data = neuron_rates_all_bins(:, baseline_bins);
        x = baseline_data(:); % Vectorize

        % Construct the post-event data (Y)
        Y = neuron_rates_all_bins(:, post_event_bins);

        % Skip if there are not enough trials or data
        if size(Y, 1) < 5 || isempty(x)
            continue;
        end

        % Prepare baseline matrix (X) for arrayROC by replicating 'x'
        X = repmat(x, 1, size(Y, 2));

        % Perform ROC analysis
        [~, ~, ~, sig] = arrayROC(X, Y, roc_perms);
        sig_results(i_neuron, :) = sig;

        % Update completed workload and cluster count
        workload_completed = workload_completed + (size(Y, 1) * ...
            size(Y, 2));
        clusters_processed = clusters_processed + 1;

        % --- Progress Reporting ---
        elapsed_time = toc(master_timer);
        rate = workload_completed / elapsed_time;
        etc_seconds = (total_workload - workload_completed) / rate;
        etc_minutes = floor(etc_seconds / 60);
        etc_rem_seconds = rem(etc_seconds, 60);

        fprintf(['Analysis completed for %d/%d clusters. ' ...
            'Approximately %d minutes and %.0f seconds remaining.\r'], ...
            clusters_processed, total_clusters_in_session, ...
            etc_minutes, etc_rem_seconds);
    end

    % Store the results for the specific condition
    analysis_results.(align_name).(cond_name).sig = sig_results;
    analysis_results.(align_name).(cond_name).time_vector = ...
        time_vector(post_event_bins);
end

if total_clusters_in_session > 0
    fprintf('\nAnalysis complete.\n');
end

end
