function selected_neurons = screen_sc_neurons(session_data, gsac_data)
% screen_sc_neurons - Identifies task-modulated neurons in the SC.
%
% This function selects neurons that show significant changes in firing rate
% across different task epochs (Baseline, Visual, Delay, Saccade). The
% logic is a re-implementation of the original local_selectNeurons.m script.
%
% INPUTS:
%   session_data - A struct containing session-specific data. It must include:
%                  - nClusters: The number of neurons.
%                  - scSide: A string ('right' or 'left') indicating the recorded SC side.
%                  - spikes.times: A cell array where each cell contains spike
%                                  timestamps for a single neuron.
%   gsac_data    - A struct containing data from the gSac task. It must include:
%                  - isHighSaliencyTrial: Logical vector for high saliency trials.
%                  - isT1Trial: Logical vector for Target 1 trials.
%                  - isT2Trial: Logical vector for Target 2 trials.
%                  - targetOnTime, fixOffTime, saccadeOnTime: Event timestamps.
%
% OUTPUT:
%   selected_neurons - A logical vector (nNeurons x 1) where true indicates
%                      a neuron that passed the selection criteria.
%

fprintf('local_selectNeurons: Identifying task-modulated neurons...\n');

% Ensure scSide is defined; default if not
if ~isfield(session_data, 'scSide')
    fprintf('WARNING in screen_sc_neurons: session_data.scSide not defined. Assuming "right" SC.\n');
    session_data.scSide = 'right';
end

% Determine contralateral target trials based on scSide
if strcmpi(session_data.scSide, 'right')      % Right SC: T1 (Left Target) is Contralateral
    selectedTrials = gsac_data.isHighSaliencyTrial & gsac_data.isT1Trial;
    fprintf('screen_sc_neurons: Right SC detected. Using T1 (Left Target) as Contralateral.\n');
elseif strcmpi(session_data.scSide, 'left') % Left SC: T2 (Right Target) is Contralateral
    selectedTrials = gsac_data.isHighSaliencyTrial & gsac_data.isT2Trial;
    fprintf('screen_sc_neurons: Left SC detected. Using T2 (Right Target) as Contralateral.\n');
else % unknown SC side
    fprintf('WARNING in screen_sc_neurons: SC side is "unknown". Using T1 as default.\n');
    selectedTrials = gsac_data.isHighSaliencyTrial & gsac_data.isT1Trial;
end

trialInd = find(selectedTrials);
nNeurons = session_data.nClusters;

% Handle edge cases: no neurons or no matching trials
if isempty(trialInd)
    fprintf('WARNING in screen_sc_neurons: No trials match selection criteria. Skipping neuron selection.\n');
    selected_neurons = false(nNeurons, 1);
    return;
end
if nNeurons == 0
    fprintf('WARNING in screen_sc_neurons: No neurons found (nClusters is 0).\n');
    selected_neurons = false(0, 1);
    return;
end

nTrials = numel(trialInd);
epoch_frs = zeros(nNeurons, 4, nTrials); % [nNeurons x nEpochs x nTrials]

% Define analysis epochs relative to key trial events
% {event, start_offset, end_offset, duration}
epochs = {
    'targetOnTime',  -0.075, 0.025,  0.1;   % 1. Baseline
    'targetOnTime',  0.05,   0.2,    0.15;  % 2. Visual
    'fixOffTime',    -0.15,  0.05,   0.2;   % 3. Delay
    'saccadeOnTime', -0.025, 0.05,   0.075  % 4. Saccade
};

% Calculate firing rates for each epoch
for i_trial = 1:nTrials
    t_idx = trialInd(i_trial);
    for i_epoch = 1:size(epochs, 1)
        event_name = epochs{i_epoch, 1};
        event_time = gsac_data.(event_name)(t_idx);

        start_offset = epochs{i_epoch, 2};
        end_offset   = epochs{i_epoch, 3};
        win_dur      = epochs{i_epoch, 4};

        if isnan(event_time)
            fprintf('Warning: NaN event time for %s on trial index %d. Setting FR to NaN.\n', event_name, t_idx);
            epoch_frs(:, i_epoch, i_trial) = NaN;
            continue;
        end

        start_time = event_time + start_offset;
        end_time   = event_time + end_offset;

        % In-line implementation of firing rate calculation (replaces local_computeFiringRate)
        for i_neuron = 1:nNeurons
            spike_times = session_data.spikes.times{i_neuron};
            spike_count = sum(spike_times >= start_time & spike_times < end_time);
            epoch_frs(i_neuron, i_epoch, i_trial) = spike_count / win_dur;
        end
    end
end

% Perform statistical tests to identify task-modulated neurons
sig_epoch_comparison = false(nNeurons, 3); % [Vis vs Base, Delay vs Base, Sac vs Base]
selected_neurons = false(nNeurons, 1);

if nTrials < 2
    fprintf('Warning: Too few trials for statistical selection. Returning empty selection.\n');
    return;
end

for i_neuron = 1:nNeurons
    % Format data for current neuron: [nTrials x nEpochs]
    neuron_frs = squeeze(epoch_frs(i_neuron, :, :))';
    % Remove trials with any NaN epochs for this neuron
    neuron_frs = neuron_frs(~any(isnan(neuron_frs), 2), :);

    if size(neuron_frs, 1) < 2
        continue; % Not enough valid trials for this neuron
    end

    try
        % Friedman test for overall modulation across epochs
        [p_friedman, ~] = friedman(neuron_frs, 1, 'off');

        if p_friedman < 0.05
            alpha_corr = 0.05 / 3; % Bonferroni correction
            comparisons = [1 2; 1 3; 1 4]; % [Base vs Vis, Base vs Delay, Base vs Sac]

            for i_comp = 1:size(comparisons, 1)
                col1 = comparisons(i_comp, 1);
                col2 = comparisons(i_comp, 2);

                % NOTE: Using Wilcoxon rank-sum test as a substitute for the missing
                % 'arrayROC' function. It's a standard non-parametric test for
                % comparing two independent samples, which serves the same
                % purpose of detecting significant differences between epochs.
                p_ranksum = ranksum(neuron_frs(:, col1), neuron_frs(:, col2));

                if p_ranksum < alpha_corr
                    sig_epoch_comparison(i_neuron, i_comp) = true;
                end
            end
        end
    catch ME
        fprintf('Statistical test failed for neuron %d: %s\n', i_neuron, ME.message);
    end

    % Select neuron if any epoch is significant AND firing rate is > 5 sp/s
    if any(sig_epoch_comparison(i_neuron, :)) && (nTrials==0 || max(mean(neuron_frs, 1, 'omitnan')) > 5)
        selected_neurons(i_neuron) = true;
    end
end

% Broaden selection criteria if too few neurons are initially selected
if nnz(selected_neurons) < min(10, nNeurons) && nNeurons > 0
    fprintf('Too few neurons selected (%d). Broadening criteria...\n', nnz(selected_neurons));
    selected_neurons = any(sig_epoch_comparison, 2);
end

end
