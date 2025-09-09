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

% Define output directory and filename
project_root = fullfile(findOneDrive, 'Code', ...
    'tokens-analysis-pipeline');
output_dir = fullfile(project_root, 'figures');

% --- Setup and Data Extraction ---
cluster_info = session_data.spikes.cluster_info;
cluster_ids = cluster_info.cluster_id;
nClusters = height(cluster_info.cluster_id);

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

% Identify gSac_jph memory-guided saccade trials with valid reward times
gSac_jph_memSac_trials = find(trialInfo.taskCode == ...
    codes.uniqueTaskCode_gSac_jph & ...
                              eventTimes.targetReillum > 0 & ...
                              eventTimes.pdsReward > 0);

% Determine which set of trials to use
if length(gSac_jph_memSac_trials) > 5 % Threshold for sufficient trials
    use_gSac_jph = true;
    memSacTrials_indices = gSac_jph_memSac_trials;
    fprintf(['Found %d valid gSac_jph memory-guided trials. ' ...
        'Using these for analysis.\n'], length(memSacTrials_indices));
else
    use_gSac_jph = false;
    % Fallback to gSac_4factors trials
    memSacTrials_logical = trialInfo.taskCode == ...
        codes.uniqueTaskCode_gSac_4factors & ...
                                ~isnan(eventTimes.targetReillum) & ...
                                eventTimes.pdsReward > 0;
    memSacTrials_indices = find(memSacTrials_logical);
    fprintf(['Insufficient gSac_jph trials. Falling back to %d ' ...
        'gSac_4factors trials.\n'], length(memSacTrials_indices));
end

if isempty(memSacTrials_indices)
    fprintf(['WARNING in screen_sc_neurons: No suitable memory ' ...
        'guided trials found. Skipping.\n']);
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

% This new neuron-centric approach iterates through each neuron and calculates
% its firing rate for all trials and all epochs in a vectorized manner,
% which is more efficient than iterating through each trial.
for i_cluster = 1:nClusters
    % Get all spike times for the current cluster
    spike_times = all_spike_times(all_spike_clusters == ...
        cluster_ids(i_cluster));

    if isempty(spike_times)
        % If no spikes, FR is 0 for all epochs and trials for this cluster
        epoch_frs(i_cluster, :, :) = 0;
        continue;
    end

    for i_epoch = 1:nEpochs
        event_name = epochs{i_epoch, 1};
        win_dur    = epochs{i_epoch, 4};

        % Get event times for all relevant trials for the current epoch
        epoch_event_times = eventTimes.(event_name)(memSacTrials_indices);

        % Create a [nMemSacTrials x 2] matrix of time windows
        time_windows = [epoch_event_times + epochs{i_epoch, 2}, ...
            epoch_event_times + epochs{i_epoch, 3}];

        % Find trials where the event time is NaN and keep track of them
        valid_trials_mask = ~isnan(epoch_event_times);

        % Filter out trials with NaN event times
        valid_time_windows = time_windows(valid_trials_mask, :);

        if isempty(valid_time_windows)
            continue; % No valid trials for this epoch, NaNs will remain
        end

        % Reshape the time_windows matrix into a single row vector of edges
        % for histcounts: [start1, end1, start2, end2, ...]
        edges = reshape(valid_time_windows', 1, []);

        try
        % Get the counts for all bins (both inside and outside the windows)
        all_counts = histcounts(spike_times, edges);
        catch me
            keyboard
        end

        % The counts within our desired windows are the odd-indexed elements
        % (1st, 3rd, 5th, etc.) of the histcounts output.
        spike_counts_in_windows = all_counts(1:2:end);

        try
        % Create a temporary array to store results for the current epoch
        temp_frs = nan(1, nMemSacTrials);
        temp_frs(valid_trials_mask) = spike_counts_in_windows / win_dur;
        catch me
            keyboard
        end

        % Place the calculated firing rates into the master matrix
        epoch_frs(i_cluster, i_epoch, :) = temp_frs;
    end
end

% --- 4. Conditional Logic for Trial Handling ---
if use_gSac_jph
    % --- Logic for gSac_jph trials ---
    fprintf(['Using gSac_jph logic to determine scSide and select ' ...
        'trials.\n']);

    % Determine scSide based on the hemifield with the most target 
    % presentations
    thetas = trialInfo.targetTheta(memSacTrials_indices)/10;
    left_vf_trials = sum(thetas > 90 & thetas < 270);
    right_vf_trials = sum(thetas < 90 | thetas > 270);

    if left_vf_trials > right_vf_trials
        scSide = 'right'; % Right SC records from the left visual field
    else
        scSide = 'left';  % Left SC records from the right visual field
    end
    fprintf(['Determined SC Side: %s based on contralateral target ' ...
        'placement.\n'], scSide);

    % Use all identified gSac_jph trials for statistical analysis
    trials_for_stats = true(1, nMemSacTrials); % Logical index for all trials

else
    % --- Logic for gSac_4factors fallback ---
    fprintf(['Using gSac_4factors logic to determine scSide and select ' ...
        'trials.\n']);

    % Determine scSide from visual epoch firing rates
    left_trials_mask = (trialInfo.targetTheta(memSacTrials_indices)/10 ...
        > 90 & trialInfo.targetTheta(memSacTrials_indices)/10 < 270);
    right_trials_mask = (trialInfo.targetTheta(memSacTrials_indices)/10 ...
        < 90 | trialInfo.targetTheta(memSacTrials_indices)/10 > 270);

    % Calculate average visual FR for left vs. right trials, for all neurons
    mean_vis_fr_left = mean(epoch_frs(:, 2, left_trials_mask), 3, ...
        'omitnan');
    mean_vis_fr_right = mean(epoch_frs(:, 2, right_trials_mask), 3, ...
        'omitnan');

    % Compare the mean across all neurons to determine side
    if mean(mean_vis_fr_left, 'omitnan') > mean(mean_vis_fr_right, ...
            'omitnan')
        scSide = 'right'; % Right SC prefers left visual field
    else
        scSide = 'left'; % Left SC prefers right visual field
    end
    fprintf(['Determined SC Side: %s by comparing visual ' ...
        'responses (L-VF vs R-VF).\n'], scSide);

    % Identify neuron-specific RFs and select trials for stats
    trials_for_stats = false(nClusters, nMemSacTrials);
    unique_locations = unique(trialInfo.targetTheta(memSacTrials_indices));

    for i_cluster = 1:nClusters
        loc_fr_sum = zeros(size(unique_locations));
        for i_loc = 1:length(unique_locations)
            loc_mask = trialInfo.targetTheta(memSacTrials_indices) == ...
                unique_locations(i_loc);
            % Sum of avg Visual FR and avg Saccade FR for this location
            vis_fr = mean(epoch_frs(i_cluster, 2, loc_mask), 3, 'omitnan');
            sac_fr = mean(epoch_frs(i_cluster, 4, loc_mask), 3, 'omitnan');
            loc_fr_sum(i_loc) = vis_fr + sac_fr;
        end
        % Find the location with the max summed FR (the RF)
        [~, rf_idx] = max(loc_fr_sum);
        rf_location(i_cluster) = unique_locations(rf_idx);
        % Mark trials at this RF location for this neuron's stats
        trials_for_stats(i_cluster, :) = (trialInfo.targetTheta( ...
            memSacTrials_indices) == rf_location(i_cluster));
    end
end

% initialize a variable to store mean firing rates across trials for each
% neuron:
all_neuron_frs = zeros(nClusters, 4);

% --- 5. Final Statistical Selection ---
for i_cluster = 1:nClusters
    if use_gSac_jph
        neuron_frs_all_trials = squeeze(epoch_frs(i_cluster, :, :))';
    else
        % Select only trials for this neuron's RF
        neuron_frs_all_trials = squeeze(epoch_frs(i_cluster, :, ...
            trials_for_stats(i_cluster, :)))';
    end

    % Remove trials with any NaN epochs
    neuron_frs = neuron_frs_all_trials(~any( ...
        isnan(neuron_frs_all_trials), 2), :);

    if size(neuron_frs, 1) < 2
        continue; % Not enough valid trials for this neuron
    else
        all_neuron_frs(i_cluster, :) = mean(neuron_frs_all_trials);
    end

    try
        [p_friedman, ~] = friedman(neuron_frs, 1, 'off');
        if p_friedman < 0.05
            alpha_corr = 0.05 / 3; % Bonferroni correction for 3 comparisons
            comparisons = [1 2; 1 3; 1 4]; % Bsl vs Vis, Bsl vs Del, Bsl vs Sac

            for i_comp = 1:size(comparisons, 1)
                p_ranksum = ranksum(neuron_frs(:, ...
                    comparisons(i_comp, 1)), neuron_frs(:, ...
                    comparisons(i_comp, 2)));
                if p_ranksum < alpha_corr
                    sig_epoch_comparison(i_cluster, i_comp) = true;
                end
            end
        end
    catch ME
        fprintf('Stat test failed for cluster %d: %s\n', ...
            cluster_ids(i_cluster), ME.message);
    end

    % Final selection criteria: significant modulation and FR > 5 Hz in 
    % any epoch
    if any(sig_epoch_comparison(i_cluster, :)) && max(mean(neuron_frs, ...
            1, 'omitnan')) > 5
        selected_neurons(i_cluster) = true;
    end
end

% generate summary figure showing mean firing rate per poch and sig epoch
% comparisons:
fig = figure('Color', 'W', 'MenuBar', 'None', 'ToolBar', 'None', ...
    'Position', ...
    [150 100 500 800]);
ax(1) = subplot(1,2,1);
imagesc(all_neuron_frs)
colormap(flipud(bone(64)))
title('Mean Firing Rate (FR)')
ax(2) = subplot(1,2,2);
imagesc(sig_epoch_comparison)
colormap(flipud(bone(64)))
title('Sig. FR Mod.')
set(ax, 'Box', 'Off', 'TickDir', 'Out')
set(ax(2), 'YTickLabel', [])
ylabel(ax(1), 'Cluster ID', 'FOntSize', 16)
set(ax(1), 'XTick', 1:4, 'XTickLabel', {'base', 'vis', 'del', 'sac'})
set(ax(2), 'XTick', 1:3, 'XTickLabel', {'vis', 'del', 'sac'})

% define PDF filename and save
figFileName = fullfile(output_dir, [session_data.metadata.unique_id, ...
    '_sc_epoch_frs.pdf']);
pdfSave(figFileName, fig.Position(3:4)/72, fig);

fprintf('Finished screening. Found %d task-modulated SC neurons.\n', nnz(selected_neurons));

end
