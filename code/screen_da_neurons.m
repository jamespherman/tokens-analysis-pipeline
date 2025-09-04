function selected_neurons = screen_da_neurons(session_data)
% screen_da_neurons - Selects putative DA neurons based on electrophysiological properties.
%
% This function screens neurons based on their baseline firing rate and
% waveform shape. It calls helper functions to calculate these metrics.
%
% INPUTS:
%   session_data - A struct containing session-specific data. It must include:
%                  - nClusters: The number of neurons.
%                  - session_id: A string identifier for the session.
%                  - spikes.wfMeans: A [nNeurons x nSamples] matrix of mean waveforms.
%                  - spikes.fs: The sampling rate of the waveforms.
%                  - (other fields required by calculate_baseline_fr)
%
% OUTPUT:
%   selected_neurons - A logical vector (nNeurons x 1) where true indicates
%                      a putative DA neuron.
%

try
fprintf('--> Screening DA neurons for session: %s\n', ...
    session_data.session_id);

nNeurons = session_data.nClusters;
if nNeurons == 0
    fprintf('WARNING: No neurons found (nClusters is 0).\n');
    selected_neurons = false(0, 1);
    return;
end

% --- 1. Calculate Metrics using Helper Functions ---

% Calculate baseline firing rate for all neurons
% This function is expected to return a vector of size [nNeurons x 1]
baseline_frs = calculate_baseline_fr(session_data);

% Calculate waveform metrics for each neuron
waveform_durations_ms = zeros(nNeurons, 1);
fs = session_data.spikes.fs;
for i = 1:nNeurons
    mean_waveform = session_data.spikes.wfMeans(i, :);
    % This function is expected to return a struct with metrics
    waveform_metrics = calculate_waveform_metrics(mean_waveform, fs);
    waveform_durations_ms(i) = waveform_metrics.peak_trough_ms;
end


% --- 2. Apply Selection Criteria ---

% Select neurons with a baseline firing rate < 20 sp/s
selected_neurons = baseline_frs < 20;

fprintf('... found %d putative DA neurons (FR < 20 sp/s).\n', ...
    nnz(selected_neurons));


% --- 3. Generate and Save Diagnostic Plot ---

% Create the 'figures' directory if it doesn't exist
if ~exist('figures', 'dir')
    mkdir('figures');
end

fig_handle = figure('Visible', 'off');
scatter(waveform_durations_ms, baseline_frs, 36, 'filled', ...
    'MarkerFaceAlpha', 0.6);
hold on;

% Highlight the selected neurons
scatter(waveform_durations_ms(selected_neurons), ...
    baseline_frs(selected_neurons), ...
        36, 'r', 'filled', 'MarkerFaceAlpha', 0.8);

% Add a line for the firing rate cutoff
yline(20, '--r', 'Label', 'FR Cutoff (20 sp/s)', 'LineWidth', 1.5);

xlabel('Waveform Peak-to-Trough Duration (ms)');
ylabel('Baseline Firing Rate (sp/s)');
title(sprintf('DA Neuron Screening: %s', session_data.session_id), ...
    'Interpreter', 'none');
legend({'All Neurons', 'Selected (Putative DA)'}, 'Location', 'northeast');
grid on;

% Save the figure
fig_filename = fullfile('figures', sprintf('da_screening_%s.png', ...
    session_data.session_id));
saveas(fig_handle, fig_filename);
fprintf('... diagnostic plot saved to %s\n', fig_filename);
close(fig_handle);

catch me
    keyboard
end

end
