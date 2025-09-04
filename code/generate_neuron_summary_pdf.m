function generate_neuron_summary_pdf(session_data, selected_neurons, unique_id)
% GENERATE_NEURON_SUMMARY_PDF Generates a multi-page PDF with diagnostic plots for each neuron.
%
%   GENERATE_NEURON_SUMMARY_PDF(session_data, selected_neurons, unique_id)
%   generates a multi-page PDF file where each page contains a summary of a
%   single neuron. The function loops through all neurons, creating a 2x2
%   panel of plots for each one.
%
%   Inputs:
%       session_data     - A struct for a single session, containing fields like
%                          'spikes' and 'cluster_info'.
%       selected_neurons - A logical vector (nNeurons x 1) indicating which
%                          neurons were selected by the screening process.
%       unique_id        - A string used to name the output PDF file.
%
%   Output:
%       The function saves a PDF file named '[unique_id]_neuron_diagnostics.pdf'
%       in the 'figures/' directory.

% Define output directory and filename
output_dir = 'figures';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
output_filename = fullfile(output_dir, [unique_id '_neuron_diagnostics.pdf']);

% Get the number of neurons
nNeurons = session_data.nClusters;

% Create a new figure for the PDF
fig = figure('Color', 'w', 'Visible', 'off', 'PaperOrientation', 'landscape');

% Pre-calculate metrics for all neurons
baseline_frs = calculate_baseline_fr(session_data);

% Loop through each neuron to create a page of plots
for i_neuron = 1:nNeurons
    % Clear the figure for the new page
    clf(fig);

    % Add a title for the page
    sgtitle(sprintf('Neuron Diagnostic Summary: Cluster %d', i_neuron));

    % --- Panel 1: Mean Waveform ---
    ax1 = mySubPlot([2, 2, 1]);
    if isfield(session_data, 'spikes') && isfield(session_data.spikes, 'wfMeans')
        plot(ax1, session_data.spikes.wfMeans(:, i_neuron));
        title(ax1, 'Mean Waveform');
        xlabel(ax1, 'Samples');
        ylabel(ax1, 'Amplitude (uV)');
        axis(ax1, 'tight');
    else
        text(0.5, 0.5, 'Waveform data not found', 'Parent', ax1, 'HorizontalAlignment', 'center');
    end

    % --- Panel 2: ISI Histogram ---
    ax2 = mySubPlot([2, 2, 2]);
    if isfield(session_data, 'spikes') && isfield(session_data.spikes, 'times')
        spike_times = session_data.spikes.times{i_neuron};
        if numel(spike_times) > 1
            isi = diff(spike_times) * 1000; % in ms
            histogram(ax2, isi, 'EdgeColor', 'k', 'FaceColor', [0.5 0.5 0.5]);
            set(ax2, 'XScale', 'log');
            title(ax2, 'ISI Histogram');
            xlabel(ax2, 'Inter-Spike Interval (ms, log scale)');
            ylabel(ax2, 'Count');
        else
            text(0.5, 0.5, 'Not enough spikes for ISI', 'Parent', ax2, 'HorizontalAlignment', 'center');
        end
    else
        text(0.5, 0.5, 'Spike times not found', 'Parent', ax2, 'HorizontalAlignment', 'center');
    end

    % --- Panel 3: Basic PSTH ---
    ax3 = mySubPlot([2, 2, 3]);
    if isfield(session_data, 'gSac') && isfield(session_data.gSac, 'outcomeOnTime')
        event_times = session_data.gSac.outcomeOnTime;
        spike_times = session_data.spikes.times{i_neuron};

        win = [-0.5, 1.0]; % Time window around event
        bin_size = 0.05; % 50 ms bins
        bins = win(1):bin_size:win(2);

        psth = zeros(size(bins));
        valid_trials = 0;
        for i_trial = 1:numel(event_times)
            if ~isnan(event_times(i_trial))
                valid_trials = valid_trials + 1;
                relative_spikes = spike_times - event_times(i_trial);
                psth = psth + histcounts(relative_spikes, [bins, bins(end)+bin_size]);
            end
        end

        if valid_trials > 0
            psth_rate = psth / (valid_trials * bin_size); % Convert to firing rate (Hz)
            bar(ax3, bins, psth_rate, 'hist');
            title(ax3, 'PSTH around Outcome');
            xlabel(ax3, 'Time from Outcome (s)');
            ylabel(ax3, 'Firing Rate (Hz)');
            xline(ax3, 0, 'r--');
        else
            text(0.5, 0.5, 'No valid trials for PSTH', 'Parent', ax3, 'HorizontalAlignment', 'center');
        end
    else
        text(0.5, 0.5, 'Event times not found', 'Parent', ax3, 'HorizontalAlignment', 'center');
    end

    % --- Panel 4: Summary Information ---
    ax4 = mySubPlot([2, 2, 4]);
    axis(ax4, 'off'); % Turn off axis for text

    % Get screening status
    if selected_neurons(i_neuron)
        screening_status = 'Selected';
    else
        screening_status = 'Not Selected';
    end

    % Get Phy quality
    if isfield(session_data, 'cluster_info') && isfield(session_data.cluster_info, 'group')
        phy_quality = session_data.cluster_info.group{i_neuron};
    else
        phy_quality = 'N/A';
    end

    % Get waveform duration
    if isfield(session_data, 'spikes') && isfield(session_data.spikes, 'fs')
        wf_metrics = calculate_waveform_metrics(session_data.spikes.wfMeans(:, i_neuron), session_data.spikes.fs);
        wf_duration = sprintf('%.2f ms', wf_metrics.peak_trough_ms);
    else
        wf_duration = 'N/A';
    end

    % Create summary text
    summary_text = {
        sprintf('Cluster ID: %d', i_neuron), ...
        sprintf('Phy Quality: %s', phy_quality), ...
        sprintf('Screening Status: %s', screening_status), ...
        sprintf('Baseline FR: %.2f Hz', baseline_frs(i_neuron)), ...
        sprintf('Waveform Duration: %s', wf_duration)
    };

    text(0.1, 0.5, summary_text, 'Parent', ax4, 'VerticalAlignment', 'middle', 'FontSize', 12);
    title(ax4, 'Summary Information');

    % Append the current figure state to the PDF
    if i_neuron == 1
        print(fig, output_filename, '-dpdf', '-fillpage');
    else
        print(fig, output_filename, '-dpdf', '-fillpage', '-append');
    end
end

% Close the figure
close(fig);

fprintf('Diagnostic PDF saved to %s\n', output_filename);

end
