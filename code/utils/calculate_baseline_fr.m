function baseline_frs = calculate_baseline_fr(session_data)
% calculate_baseline_fr - Computes the average baseline firing rate for neurons.
%
% This function calculates the firing rate for each neuron during a
% predefined baseline period, averaged across all valid trials.
%
% INPUTS:
%   session_data - A struct containing session-specific data, conforming to
%                  the structure defined in 'docs/preprocessing_docs/session_data_dictionary.md'.
%                  It must include:
%                  - spikes.cluster_info: Table with information about each cluster.
%                  - spikes.times: Vector of all spike times.
%                  - spikes.clusters: Vector mapping each spike to a cluster ID.
%                  - eventTimes.targOn: Vector of target onset times for each trial.
%
% OUTPUT:
%   baseline_frs - A vector (nClusters x 1) containing the average baseline
%                  firing rate for each cluster in spikes/sec.
%

% --- Setup ---
% Get cluster information from the standardized session_data structure.
if ~isfield(session_data, 'spikes') || ~isfield(session_data.spikes, 'cluster_info')
    error('calculate_baseline_fr: session_data.spikes.cluster_info not found.');
end
cluster_info = session_data.spikes.cluster_info;
cluster_ids = cluster_info.cluster_id;
nClusters = height(cluster_info);

% Get all spike times and their cluster assignments.
all_spike_times = session_data.spikes.times;
all_spike_clusters = session_data.spikes.clusters;

% Get event times from the standardized 'eventTimes' struct.
if isfield(session_data, 'eventTimes') && isfield(session_data.eventTimes, 'targOn')
    target_on_times = session_data.eventTimes.targOn;
else
    error('calculate_baseline_fr: session_data.eventTimes.targOn not found.');
end

% --- Define Baseline Period ---
% We define the baseline period as the 100ms window immediately preceding
% the target onset on each trial.
baseline_start_offset = -0.100; % 100 ms before target onset
baseline_end_offset   = 0.000;  % at target onset
baseline_duration_sec = baseline_end_offset - baseline_start_offset; % should be 0.1

% --- Calculate Firing Rates ---
baseline_frs = zeros(nClusters, 1);
valid_trials = find(~isnan(target_on_times));
nValidTrials = numel(valid_trials);

if nValidTrials == 0
    fprintf('WARNING in calculate_baseline_fr: No valid trials with target onset found.\n');
    return; % Returns a vector of zeros
end

% Total duration across all valid trials, used for averaging.
total_baseline_duration = nValidTrials * baseline_duration_sec;

if total_baseline_duration <= 0
    fprintf('WARNING in calculate_baseline_fr: Total baseline duration is zero or negative.\n');
    return; % Returns a vector of zeros
end

% Iterate through each cluster and calculate its baseline firing rate.
for i_cluster = 1:nClusters
    cluster_id = cluster_ids(i_cluster);

    % Get spike times for the current cluster.
    neuron_spike_times = all_spike_times(all_spike_clusters == cluster_id);

    if isempty(neuron_spike_times)
        baseline_frs(i_cluster) = 0;
        continue;
    end

    total_spike_count = 0;

    % Sum spikes across the baseline window of all valid trials.
    for i_trial = 1:nValidTrials
        trial_idx = valid_trials(i_trial);
        event_time = target_on_times(trial_idx);

        start_time = event_time + baseline_start_offset;
        end_time   = event_time + baseline_end_offset;

        % Count spikes within the baseline window for this trial.
        total_spike_count = total_spike_count + sum(neuron_spike_times >= start_time & neuron_spike_times < end_time);
    end

    % Calculate average firing rate for this neuron.
    baseline_frs(i_cluster) = total_spike_count / total_baseline_duration;
end

end
