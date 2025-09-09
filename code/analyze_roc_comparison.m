%% analyze_roc_comparison.m
%
%   Performs a bin-by-bin ROC analysis to compare firing rates between two
%   specified conditions for key trial events.
%
% INPUTS:
%   core_data  - A struct containing processed neuronal data.
%   conditions - A struct with logical masks for different trial conditions.
%
% OUTPUT:
%   analysis_results - A struct containing the significance matrix for each
%                      of the three main comparisons.
%
% Author: Jules
% Date: 2025-09-09

function analysis_results = analyze_roc_comparison(core_data, conditions)

%% Setup Paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Analysis Parameters
roc_perms = 200;
n_neurons = size(core_data.spikes.CUE_ON.rates, 1);

% Define the three key comparisons to make
comparisons = {
    ... % Comparison 1: Normal vs. Uniform at Cue On
    {'CUE_ON', 'is_normal_dist', 'is_uniform_dist', 'cue_normal_vs_uniform'}, ...
    ... % Comparison 2: Common vs. Rare High at Outcome
    {'outcomeOn', 'is_norm_common', 'is_norm_rare_high', 'outcome_common_vs_rare'}, ...
    ... % Comparison 3: Common vs. Rare High at Reward
    {'reward', 'is_common_reward_no_spe', 'is_rare_high_reward_no_spe', 'reward_common_vs_rare'}
    };

%% Main Analysis Loop
% Iterate through each of the three comparisons
for i_comp = 1:numel(comparisons)
    comp_params = comparisons{i_comp};
    align_name = comp_params{1};
    cond1_name = comp_params{2};
    cond2_name = comp_params{3};
    result_name = comp_params{4};

    % Get data for the current alignment
    rates_all_neurons = core_data.spikes.(align_name).rates;
    time_vector = core_data.spikes.(align_name).time_vector;
    n_time_bins = numel(time_vector);

    % Get trial indices for the two conditions
    trials_cond1 = conditions.(cond1_name);
    trials_cond2 = conditions.(cond2_name);

    % Preallocate storage for significance results
    sig_results = nan(n_neurons, n_time_bins);

    % Iterate through each neuron
    for i_neuron = 1:n_neurons
        % Extract rates for the current neuron, all trials, all bins
        neuron_rates_all_trials = squeeze(rates_all_neurons(i_neuron, :, :));

        % Get the trial-by-bin matrices for each condition
        rates_cond1 = neuron_rates_all_trials(trials_cond1, :);
        rates_cond2 = neuron_rates_all_trials(trials_cond2, :);

        % Skip if there are not enough trials in either condition
        if size(rates_cond1, 1) < 5 || size(rates_cond2, 1) < 5
            continue;
        end

        % Perform bin-by-bin ROC analysis
        % Note: arrayROC is efficient for this; no need to loop bins
        [~, ~, ~, sig] = arrayROC(rates_cond1, rates_cond2, roc_perms);
        sig_results(i_neuron, :) = sig;
    end

    % Store the results
    analysis_results.(result_name).sig = sig_results;
    analysis_results.(result_name).time_vector = time_vector;
    analysis_results.(result_name).cond_names = {cond1_name, cond2_name};
end

end
