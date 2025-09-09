%% analyze_baseline_comparison.m
%
%   Compares baseline firing rates to post-event firing rates for specified
%   conditions and alignment events using ROC analysis.
%
% INPUTS:
%   core_data  - A struct containing processed neuronal data.
%   conditions - A struct with logical masks for different trial conditions.
%
% OUTPUT:
%   analysis_results - A struct containing the significance matrix for each
%                      condition and alignment.
%
% Author: Jules
% Date: 2025-09-09

function analysis_results = analyze_baseline_comparison(core_data, conditions)

%% Setup Paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Analysis Parameters
baseline_window = [-0.2, 0]; % -200ms to 0ms relative to event
roc_perms = 200;

% Define the alignment events and conditions to be analyzed
align_events = {'CUE_ON', 'outcomeOn', 'reward'};
conds_to_analyze = {
    'is_normal_dist', ...
    'is_uniform_dist', ...
    'is_norm_common', ...
    'is_norm_rare_high', ...
    'is_common_reward_no_spe', ...
    'is_rare_high_reward_no_spe'
    };

n_neurons = size(core_data.spikes.CUE_ON.rates, 1);

%% Main Analysis Loop
% Iterate through each alignment event
for i_align = 1:numel(align_events)
    align_name = align_events{i_align};
    time_vector = core_data.spikes.(align_name).time_vector;

    % Find time bins corresponding to the baseline and post-event periods
    baseline_bins = time_vector >= baseline_window(1) & ...
        time_vector < baseline_window(2);
    post_event_bins = time_vector >= 0;

    % Iterate through each condition
    for i_cond = 1:numel(conds_to_analyze)
        cond_name = conds_to_analyze{i_cond};
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
        end

        % Store the results
        analysis_results.(align_name).(cond_name).sig = sig_results;
        analysis_results.(align_name).time_vector = ...
            time_vector(post_event_bins);
    end
end

end
