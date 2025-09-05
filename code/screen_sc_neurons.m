function [selected_neurons, sig_epoch_comparison, scSide] = screen_sc_neurons(session_data)
% screen_sc_neurons - Identifies task-modulated neurons in the SC.
%
% This function selects neurons that show significant changes in firing rate
% across different task epochs (Baseline, Visual, Delay, Saccade). It uses
% the standardized `session_data` structure. It automatically determines
% the recorded side of the SC based on neural activity.
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

% --- Setup and Data Extraction ---
cluster_info = session_data.spikes.cluster_info;
cluster_ids = cluster_info.cluster_id;
nClusters = height(cluster_info);

all_spike_times = session_data.spikes.times;
all_spike_clusters = session_data.spikes.clusters;

trialInfo = session_data.trialInfo;
eventTimes = session_data.eventTimes;

% Initialize outputs for early return
selected_neurons = false(nClusters, 1);
sig_epoch_comparison = false(nClusters, 3);
scSide = 'unknown';

if nClusters == 0
    fprintf('WARNING in screen_sc_neurons: No clusters found.\n');
    return;
end

% --- Trial Selection ---
% Reformat isVisSac from cell to logical
tempIsVisSac = false(height(trialInfo), 1);
for i = 1:length(trialInfo.isVisSac)
    if ~isempty(trialInfo.isVisSac{i})
        tempIsVisSac(i) = logical(trialInfo.isVisSac{i});
    end
end
trialInfo.isVisSac = tempIsVisSac;

% Identify memory guided saccade trials from gSac_jph task
codes = initCodes();
gSac_jph_trials = trialInfo.taskCode == codes.uniqueTaskCode_gSac_jph;
memSacTrials = gSac_jph_trials & ~trialInfo.isVisSac;

% Fallback if no memory guided trials in gSac_jph
if nnz(memSacTrials) == 0
    fprintf('No memory guided gSac_jph trials found. Using all rewarded memory guided trials.\n');
    memSacTrials = ~trialInfo.isVisSac & eventTimes.reward > 0 & eventTimes.targetReillum > 0;
end

if nnz(memSacTrials) == 0
    fprintf('WARNING in screen_sc_neurons: No memory guided trials found. Skipping.\n');
    return;
end

% --- Determine SC Side from Neural Activity ---
left_trials = find(memSacTrials & trialInfo.targTheta_x10 > 0);
right_trials = find(memSacTrials & trialInfo.targTheta_x10 < 0);

if isempty(left_trials) || isempty(right_trials)
    fprintf('WARNING in screen_sc_neurons: Not enough left/right trials to determine SC side.\n');
    return;
end

% --- Epoch Definition ---
epochs = {
    'targetOn', -0.075, 0.025, 0.1;   % 1. Baseline
    'targetOn', 0.05,   0.2,   0.15;  % 2. Visual
    'fixOff',   -0.15,  0.05,  0.2;   % 3. Delay
    'saccadeOnset', -0.025, 0.05,  0.075  % 4. Saccade
};

% --- Firing Rate Calculation for SC Side Determination ---
% Use a simplified visual epoch for side determination
vis_epoch_fr = zeros(nClusters, nnz(memSacTrials));
all_mem_trials = find(memSacTrials);

for i_trial = 1:length(all_mem_trials)
    t_idx = all_mem_trials(i_trial);
    event_time = eventTimes.targetOn(t_idx);
    if isnan(event_time)
        vis_epoch_fr(:, i_trial) = NaN;
        continue;
    end
    start_time = event_time + 0.05; % Visual epoch start
    end_time = event_time + 0.2;   % Visual epoch end
    win_dur = 0.15;

    for i_cluster = 1:nClusters
        cluster_id = cluster_ids(i_cluster);
        spike_times = all_spike_times(all_spike_clusters == cluster_id);
        spike_count = sum(spike_times >= start_time & spike_times < end_time);
        vis_epoch_fr(i_cluster, i_trial) = spike_count / win_dur;
    end
end

left_fr = mean(vis_epoch_fr(:, ismember(all_mem_trials, left_trials)), 2, 'omitnan');
right_fr = mean(vis_epoch_fr(:, ismember(all_mem_trials, right_trials)), 2, 'omitnan');

if mean(left_fr) > mean(right_fr)
    scSide = 'right'; % Right SC prefers left visual field
    is_contralateral = trialInfo.targTheta_x10 > 0;
    fprintf('screen_sc_neurons: Determined SC Side: right (L-VF > R-VF)\n');
else
    scSide = 'left'; % Left SC prefers right visual field
    is_contralateral = trialInfo.targTheta_x10 < 0;
    fprintf('screen_sc_neurons: Determined SC Side: left (R-VF > L-VF)\n');
end

trialInd = find(memSacTrials & is_contralateral);
nTrials = numel(trialInd);

if nTrials < 2
    fprintf('Warning: Too few contralateral trials (%d) for statistical selection.\n', nTrials);
    return;
end

epoch_frs = zeros(nClusters, 4, nTrials); % [nClusters x nEpochs x nTrials]

% --- Firing Rate Calculation for Selected Trials ---
for i_trial = 1:nTrials
    t_idx = trialInd(i_trial);
    for i_epoch = 1:size(epochs, 1)
        event_name = epochs{i_epoch, 1};
        event_time = eventTimes.(event_name)(t_idx);

        start_offset = epochs{i_epoch, 2};
        end_offset   = epochs{i_epoch, 3};
        win_dur      = epochs{i_epoch, 4};

        if isnan(event_time)
            epoch_frs(:, i_epoch, i_trial) = NaN;
            continue;
        end

        start_time = event_time + start_offset;
        end_time   = event_time + end_offset;

        for i_cluster = 1:nClusters
            cluster_id = cluster_ids(i_cluster);
            spike_times = all_spike_times(all_spike_clusters == cluster_id);
            spike_count = sum(spike_times >= start_time & spike_times < end_time);
            epoch_frs(i_cluster, i_epoch, i_trial) = spike_count / win_dur;
        end
    end
end

% --- Statistical Selection ---
for i_cluster = 1:nClusters
    neuron_frs = squeeze(epoch_frs(i_cluster, :, :))';
    neuron_frs = neuron_frs(~any(isnan(neuron_frs), 2), :);

    if size(neuron_frs, 1) < 2
        continue; % Not enough valid trials for this neuron
    end

    try
        [p_friedman, ~] = friedman(neuron_frs, 1, 'off');
        if p_friedman < 0.05
            alpha_corr = 0.05 / 3; % Bonferroni correction
            comparisons = [1 2; 1 3; 1 4];

            for i_comp = 1:size(comparisons, 1)
                col1 = comparisons(i_comp, 1);
                col2 = comparisons(i_comp, 2);
                p_ranksum = ranksum(neuron_frs(:, col1), neuron_frs(:, col2));
                if p_ranksum < alpha_corr
                    sig_epoch_comparison(i_cluster, i_comp) = true;
                end
            end
        end
    catch ME
        fprintf('Stat test failed for cluster %d: %s\n', cluster_ids(i_cluster), ME.message);
    end

    if any(sig_epoch_comparison(i_cluster, :)) && max(mean(neuron_frs, 1, 'omitnan')) > 5
        selected_neurons(i_cluster) = true;
    end
end

% Broaden selection criteria if too few neurons are initially selected
if nnz(selected_neurons) < min(10, nClusters) && nClusters > 0
    fprintf('Too few neurons selected (%d). Broadening criteria...\n', nnz(selected_neurons));
    selected_neurons = any(sig_epoch_comparison, 2);
end

end
