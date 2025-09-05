function [selected_neurons, sig_epoch_comparison, scSide] = screen_sc_neurons(session_data)
% screen_sc_neurons - Identifies task-modulated neurons in the SC using a
% neuron-centric, vectorized approach.
%
% This function selects neurons that show significant changes in firing rate
% across different task epochs (Baseline, Visual, Delay, Saccade). It first
% attempts to use memory-guided saccade trials from the 'gSac_jph' task.
% If insufficient trials are available, it falls back to using trials from
% the 'gSac_4factors' task. The logic for determining the recorded side of
% the SC and for selecting trials for analysis differs based on the task used.
%
% INPUTS:
%   session_data - A struct containing session-specific data, conforming to the
%                  `session_data_dictionary.md`.
%
% OUTPUT:
%   selected_neurons - A logical vector (nClusters x 1) where true indicates
%                      a neuron that passed the selection criteria.
%   sig_epoch_comparison - A logical matrix (nClusters x 3) indicating
%                          significant firing rate changes between epochs.
%   scSide           - A string ('right' or 'left') indicating the determined
%                      recorded SC side.
%

fprintf('screen_sc_neurons: Identifying task-modulated neurons...\n');

% --- 1. Setup and Data Extraction ---
cluster_info = session_data.spikes.cluster_info;
cluster_ids = cluster_info.cluster_id;
nClusters = height(cluster_info);

all_spike_times = session_data.spikes.times;
all_spike_clusters = session_data.spikes.clusters;

trialInfo = session_data.trialInfo;
eventTimes = session_data.eventTimes;
codes = initCodes();

% Initialize outputs for early return
selected_neurons = false(nClusters, 1);
sig_epoch_comparison = false(nClusters, 3);
scSide = 'unknown';

if nClusters == 0
    fprintf('WARNING in screen_sc_neurons: No clusters found.\n');
    return;
end

% --- 2. Initial Trial Identification ---
% Reformat isVisSac from cell to logical if necessary
if iscell(trialInfo.isVisSac)
    tempIsVisSac = false(height(trialInfo), 1);
    for i = 1:length(trialInfo.isVisSac)
        if ~isempty(trialInfo.isVisSac{i})
            tempIsVisSac(i) = logical(trialInfo.isVisSac{i});
        end
    end
    trialInfo.isVisSac = tempIsVisSac;
end

% Identify all memory-guided saccade trials
all_memSac_trials = find(~trialInfo.isVisSac);

% Identify gSac_jph memory-guided saccade trials with valid reward times
gSac_jph_memSac_trials = find(trialInfo.taskCode == codes.uniqueTaskCode_gSac_jph & ...
                              ~trialInfo.isVisSac & ...
                              ~isnan(eventTimes.reward));

% Determine which set of trials to use
if length(gSac_jph_memSac_trials) > 5 % Threshold for sufficient trials
    use_gSac_jph = true;
    memSacTrials_indices = gSac_jph_memSac_trials;
    fprintf('Found %d valid gSac_jph memory-guided trials. Using these for analysis.\n', length(memSacTrials_indices));
else
    use_gSac_jph = false;
    % Fallback to gSac_4factors trials
    memSacTrials_indices = find(trialInfo.taskCode == codes.uniqueTaskCode_gSac_4factors & ...
                                ~trialInfo.isVisSac & ...
                                ~isnan(eventTimes.reward));
    fprintf('Insufficient gSac_jph trials. Falling back to %d gSac_4factors trials.\n', length(memSacTrials_indices));
end

if isempty(memSacTrials_indices)
    fprintf('WARNING in screen_sc_neurons: No suitable memory guided trials found. Skipping.\n');
    return;
end

nMemSacTrials = length(memSacTrials_indices);

% --- 3. Vectorized Firing Rate Calculation ---
% Epoch definitions: {event_name, start_offset, end_offset, duration}
epochs = {
    'targetOn',     -0.075, 0.025,  0.1;   % 1. Baseline
    'targetOn',     0.05,   0.2,    0.15;  % 2. Visual
    'fixOff',       -0.15,  0.05,   0.2;   % 3. Delay
    'saccadeOnset', -0.025, 0.05,   0.075  % 4. Saccade
};
nEpochs = size(epochs, 1);
epoch_frs = nan(nClusters, nEpochs, nMemSacTrials);

% Bin edges for histcounts, adding an extra edge for the last cluster
cluster_id_edges = [cluster_ids; cluster_ids(end)+1];

for i_trial = 1:nMemSacTrials
    t_idx = memSacTrials_indices(i_trial);
    for i_epoch = 1:nEpochs
        event_name = epochs{i_epoch, 1};
        event_time = eventTimes.(event_name)(t_idx);

        if isnan(event_time)
            continue; % epoch_frs remains NaN
        end

        start_time = event_time + epochs{i_epoch, 2};
        end_time   = event_time + epochs{i_epoch, 3};
        win_dur    = epochs{i_epoch, 4};

        % Find spikes within the epoch window
        spikes_in_epoch = all_spike_times >= start_time & all_spike_times < end_time;
        clusters_in_epoch = all_spike_clusters(spikes_in_epoch);

        % Count spikes per cluster in a vectorized way
        spike_counts = histcounts(clusters_in_epoch, cluster_id_edges);

        % Calculate and store firing rate
        epoch_frs(:, i_epoch, i_trial) = spike_counts / win_dur;
    end
end

% --- 4. Conditional Logic for Trial Handling ---
if use_gSac_jph
    % --- Logic for gSac_jph trials ---
    fprintf('Using gSac_jph logic to determine scSide and select trials.\n');

    % Determine scSide directly from target locations
    avg_target_theta = mean(trialInfo.targTheta_x10(memSacTrials_indices), 'omitnan');
    if avg_target_theta > 0 % Left visual field
        scSide = 'right';
    else % Right visual field
        scSide = 'left';
    end
    fprintf('Determined SC Side: %s based on contralateral target placement.\n', scSide);

    % Use all identified gSac_jph trials for statistical analysis
    trials_for_stats = true(1, nMemSacTrials); % Logical index for all trials

else
    % --- Logic for gSac_4factors fallback ---
    fprintf('Using gSac_4factors logic to determine scSide and select trials.\n');

    % Determine scSide from visual epoch firing rates
    left_trials_mask = trialInfo.targTheta_x10(memSacTrials_indices) > 0;
    right_trials_mask = trialInfo.targTheta_x10(memSacTrials_indices) < 0;

    % Calculate average visual FR for left vs. right trials, for all neurons
    mean_vis_fr_left = mean(epoch_frs(:, 2, left_trials_mask), 3, 'omitnan');
    mean_vis_fr_right = mean(epoch_frs(:, 2, right_trials_mask), 3, 'omitnan');

    % Compare the mean across all neurons to determine side
    if mean(mean_vis_fr_left, 'omitnan') > mean(mean_vis_fr_right, 'omitnan')
        scSide = 'right'; % Right SC prefers left visual field
    else
        scSide = 'left'; % Left SC prefers right visual field
    end
    fprintf('Determined SC Side: %s by comparing visual responses (L-VF vs R-VF).\n', scSide);

    % Identify neuron-specific RFs and select trials for stats
    trials_for_stats = false(nClusters, nMemSacTrials);
    unique_locations = unique(trialInfo.targTheta_x10(memSacTrials_indices));

    for i_cluster = 1:nClusters
        loc_fr_sum = zeros(size(unique_locations));
        for i_loc = 1:length(unique_locations)
            loc_mask = trialInfo.targTheta_x10(memSacTrials_indices) == unique_locations(i_loc);
            % Sum of avg Visual FR and avg Saccade FR for this location
            vis_fr = mean(epoch_frs(i_cluster, 2, loc_mask), 3, 'omitnan');
            sac_fr = mean(epoch_frs(i_cluster, 4, loc_mask), 3, 'omitnan');
            loc_fr_sum(i_loc) = vis_fr + sac_fr;
        end
        % Find the location with the max summed FR (the RF)
        [~, rf_idx] = max(loc_fr_sum);
        rf_location = unique_locations(rf_idx);
        % Mark trials at this RF location for this neuron's stats
        trials_for_stats(i_cluster, :) = (trialInfo.targTheta_x10(memSacTrials_indices) == rf_location);
    end
end

% --- 5. Final Statistical Selection ---
for i_cluster = 1:nClusters
    if use_gSac_jph
        neuron_frs_all_trials = squeeze(epoch_frs(i_cluster, :, :))';
    else
        % Select only trials for this neuron's RF
        neuron_frs_all_trials = squeeze(epoch_frs(i_cluster, :, trials_for_stats(i_cluster, :)))';
    end

    % Remove trials with any NaN epochs
    neuron_frs = neuron_frs_all_trials(~any(isnan(neuron_frs_all_trials), 2), :);

    if size(neuron_frs, 1) < 2
        continue; % Not enough valid trials for this neuron
    end

    try
        [p_friedman, ~] = friedman(neuron_frs, 1, 'off');
        if p_friedman < 0.05
            alpha_corr = 0.05 / 3; % Bonferroni correction for 3 comparisons
            comparisons = [1 2; 1 3; 1 4]; % Bsl vs Vis, Bsl vs Del, Bsl vs Sac

            for i_comp = 1:size(comparisons, 1)
                p_ranksum = ranksum(neuron_frs(:, comparisons(i_comp, 1)), neuron_frs(:, comparisons(i_comp, 2)));
                if p_ranksum < alpha_corr
                    sig_epoch_comparison(i_cluster, i_comp) = true;
                end
            end
        end
    catch ME
        fprintf('Stat test failed for cluster %d: %s\n', cluster_ids(i_cluster), ME.message);
    end

    % Final selection criteria: significant modulation and FR > 5 Hz in any epoch
    if any(sig_epoch_comparison(i_cluster, :)) && max(mean(neuron_frs, 1, 'omitnan')) > 5
        selected_neurons(i_cluster) = true;
    end
end

fprintf('Finished screening. Found %d task-modulated SC neurons.\n', nnz(selected_neurons));

end
