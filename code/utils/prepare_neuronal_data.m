%% prepare_neuronal_data.m
%
%   Generates binned, event-aligned firing rate matrices for selected
%   neurons.
%
% INPUTS:
%   session_data         - The main data structure.
%   selected_neurons     - A logical vector for neuron selection.
%   tokens_trial_indices - Indices for trials of the 'tokens' task.
%
% OUTPUT:
%   aligned_spikes       - A struct with aligned spike data.
%
% Author: Jules
% Date: 2025-09-08

function aligned_spikes = prepare_neuronal_data(session_data, ...
    selected_neurons, tokens_trial_indices)

%% Define Alignment Parameters
alignment_events = {'cueOn', 'outcomeOn'};
time_window = [-0.5, 1.5]; % in seconds
bin_width = 0.025; % 25ms in seconds

% Get the list of neuron cluster IDs to analyze
neuron_cluster_ids = find(selected_neurons);
n_selected_neurons = numel(neuron_cluster_ids);
n_tokens_trials = numel(tokens_trial_indices);

%% Process Each Alignment Event
for i_event = 1:numel(alignment_events)
    event_name = alignment_events{i_event};

    % Get the event times for the current alignment event, for tokens trials
    % only
    event_times = session_data.behavior.eventTimes.(event_name)...
        (tokens_trial_indices);

    % Initialize storage for the binned spike rates for this event
    [~, ~, time_vector] = alignAndBinSpikes([], [], time_window(1), ...
        time_window(2), bin_width);
    n_time_bins = numel(time_vector);

    binned_rates = nan(n_selected_neurons, n_tokens_trials, n_time_bins);

    % Loop through each selected neuron
    for i_neuron = 1:n_selected_neurons
        cluster_id = neuron_cluster_ids(i_neuron);

        % Get spike times for the current neuron
        spike_times = session_data.spikes.times(session_data.spikes.clusters == cluster_id);

        % Align and bin spikes for the current neuron and event
        [~, binned_counts, ~] = alignAndBinSpikes(spike_times, ...
            event_times, time_window(1), time_window(2), bin_width);

        % Convert counts to firing rate (spikes/sec)
        binned_rates(i_neuron, :, :) = binned_counts / bin_width;
    end

    % Store the results for this event
    aligned_spikes.(event_name).rates = binned_rates;
    aligned_spikes.(event_name).time_vector = time_vector;
end

end
