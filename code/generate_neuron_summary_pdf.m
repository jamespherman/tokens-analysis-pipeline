function generate_neuron_summary_pdf(session_data, selected_neurons, unique_id)
% GENERATE_NEURON_SUMMARY_PDF Generates a multi-page PDF with diagnostic plots for each neuron.
%
%   GENERATE_NEURON_SUMMARY_PDF(session_data, selected_neurons, unique_id)
%   generates a multi-page PDF file where each page contains a summary of a
%   single neuron, conforming to the standardized `session_data` structure.
%
%   Inputs:
%       session_data     - A struct for a single session, containing fields like
%                          'spikes' and 'cluster_info'.
%       selected_neurons - A logical vector (nClusters x 1) indicating which
%                          neurons were selected by a screening process.
%       unique_id        - A string used to name the output PDF file.
%
%   Output:
%       The function saves a PDF file named '[unique_id]_neuron_diagnostics.pdf'
%       in the 'figures/' directory.

% Define output directory and filename
project_root = fullfile(findOneDrive, 'Code', 'tokens-analysis-pipeline');
output_dir = fullfile(project_root, 'figures');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
output_filename = fullfile(output_dir, [unique_id '_neuron_diagnostics.pdf']);

% --- Setup and Data Extraction ---
cluster_info = session_data.spikes.cluster_info;
cluster_ids = cluster_info.cluster_id;
nClusters = height(cluster_info);
all_spike_times = session_data.spikes.times;
all_spike_clusters = session_data.spikes.clusters;

% Create a new figure for the PDF
fig = figure('Color', 'w', 'Visible', 'off', 'PaperOrientation', 'landscape');

% Loop through each cluster to create a page of plots
for i_cluster = 1:nClusters
    cluster_id = cluster_ids(i_cluster);

    % Clear the figure for the new page
    clf(fig);

    % Add a title for the page
    sgtitle(sprintf('Neuron Diagnostic Summary: Cluster %d', cluster_id));

    % --- Panel 1: Mean Waveform ---
    ax1 = mySubPlot([2, 2, 1]);
    if isfield(session_data.spikes, 'wfMeans') && numel(session_data.spikes.wfMeans) >= i_cluster
        plot(ax1, session_data.spikes.wfMeans{i_cluster}');
        title(ax1, 'Mean Waveform');
        xlabel(ax1, 'Samples');
        ylabel(ax1, 'Amplitude (uV)');
        axis(ax1, 'tight');
    else
        text(0.5, 0.5, 'Waveform data not found', 'Parent', ax1, 'HorizontalAlignment', 'center');
    end

    % --- Panel 2: ISI Histogram ---
    ax2 = mySubPlot([2, 2, 2]);
    spike_times = all_spike_times(all_spike_clusters == cluster_id);
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

    % --- Panel 3: Basic PSTH ---
    ax3 = mySubPlot([2, 2, 3]);
    if isfield(session_data, 'eventTimes') && isfield(session_data.eventTimes, 'pdsOutcomeOn')
        event_times = session_data.eventTimes.pdsOutcomeOn;
        win = [-0.5, 1.0]; % Time window around event
        bin_size = 0.05; % 50 ms bins
        bins = win(1):bin_size:win(2);

        valid_trials_idx = find(~isnan(event_times));
        n_valid_trials = numel(valid_trials_idx);
        psth = zeros(size(bins));

        if n_valid_trials > 0
            for i_trial = 1:n_valid_trials
                event_time = event_times(valid_trials_idx(i_trial));
                relative_spikes = spike_times - event_time;
                psth = psth + histcounts(relative_spikes, [bins, bins(end)+bin_size]);
            end
            psth_rate = psth / (n_valid_trials * bin_size); % Convert to firing rate (Hz)
            bar(ax3, bins, psth_rate, 'hist');
            title(ax3, 'PSTH around Outcome');
            xlabel(ax3, 'Time from Outcome (s)');
            ylabel(ax3, 'Firing Rate (Hz)');
            xline(ax3, 0, 'r--');
        else
            text(0.5, 0.5, 'No valid trials for PSTH', 'Parent', ...
                ax3, 'HorizontalAlignment', 'center');
        end
    else
        text(0.5, 0.5, 'Event times (pdsOutcomeOn) not found', ...
            'Parent', ax3, 'HorizontalAlignment', 'center');
    end

    % --- Panel 4: Summary Information ---
    ax4 = mySubPlot([2, 2, 4]);
    axis(ax4, 'off'); % Turn off axis for text

    screening_status = 'Not Selected';
    if numel(selected_neurons) >= i_cluster && selected_neurons(i_cluster)
        screening_status = 'Selected';
    end

    phy_quality = 'N/A';
    % Find the row corresponding to the current cluster_id
    info_row = cluster_info.cluster_id == cluster_id;
    if any(info_row) && isfield(cluster_info, 'group')
        phy_quality = cluster_info.group(info_row);
    end

    % Retrieve pre-calculated metrics
    baseline_fr = session_data.metrics.baseline_frs(i_cluster);
    wf_metrics = session_data.metrics.wf_metrics(i_cluster);
    wf_duration = sprintf('%.2f ms', wf_metrics.peak_trough_ms);

    % Create summary text
    summary_text = {
        sprintf('Cluster ID: %d', cluster_id), ...
        sprintf('Phy Quality: %s', phy_quality), ...
        sprintf('Screening Status: %s', screening_status), ...
        sprintf('Baseline FR: %.2f Hz', baseline_fr), ...
        sprintf('Waveform Duration: %s', wf_duration)
    };

    text(0.1, 0.5, summary_text, 'Parent', ax4, ...
        'VerticalAlignment', 'middle', 'FontSize', 12);
    title(ax4, 'Summary Information');

    % Append the current figure state to the PDF
    if i_cluster == 1
        print(fig, output_filename, '-dpdf', '-fillpage');
    else
        print(fig, output_filename, '-dpdf', '-fillpage', '-append');
    end
end

% Close the figure
close(fig);

fprintf('Diagnostic PDF saved to %s\n', output_filename);

end