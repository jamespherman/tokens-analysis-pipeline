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

% In-line function to report timing
giveFeed = @(x)disp([num2str(toc) 's - ' x]);

% Define output directory and filename
project_root = fullfile(findOneDrive, 'Code', 'tokens-analysis-pipeline');
output_dir = fullfile(project_root, 'figures');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end


% --- Setup and Data Extraction ---
cluster_info = session_data.spikes.cluster_info;
cluster_ids = cluster_info.cluster_id;
nClusters = height(cluster_info.cluster_id);
all_spike_times = session_data.spikes.times;
all_spike_clusters = session_data.spikes.clusters;

% Create a new figure for the PDF
fig = figure('Color', 'w', 'Visible', 'off', 'PaperOrientation', ...
    'landscape');

% --- Setup and Data Extraction ---
codes = initCodes;
gSac_task_code = codes.uniqueTaskCode_gSac_4factors;
has_gSac_trials = any(session_data.trialInfo.taskCode == gSac_task_code);
has_gSac_events = isfield(session_data.eventTimes, 'targetOn') && ...
                  isfield(session_data.eventTimes, 'saccadeOnset');
do_spatial_tuning_plot = has_gSac_trials && has_gSac_events;


% Loop through each cluster to create a page of plots
for i_cluster = 1:nClusters

    % define CLUSTER ID:
    cluster_id = cluster_ids(i_cluster);
    
    % define output file name:
    output_filename = fullfile(output_dir, ...
        [unique_id '_' sprintf('%03d',i_cluster) ...
        '_neuron_diagnostics.pdf']);

    % Clear the figure for the new page
    clf(fig);

    % Add a title for the page
    sgtitle(sprintf('Neuron Diagnostic Summary: Cluster %d', cluster_id));

    spike_times = all_spike_times(all_spike_clusters == cluster_id);

    if do_spatial_tuning_plot
        % --- Create the 12-panel plot for gSac_4factors sessions ---

        % --- Top Row: Original Diagnostics ---
        % Panel 1: Mean Waveform
        ax1 = subplot(6, 6, [1, 2]);
        if isfield(session_data.spikes, 'wfMeans') && ...
            numel(session_data.spikes.wfMeans) >= i_cluster
            plot(ax1, session_data.spikes.wfMeans{i_cluster}');
            title(ax1, 'Mean Waveform');
            xlabel(ax1, 'Samples');
            ylabel(ax1, 'Amplitude (uV)');
            axis(ax1, 'tight');
        else
            text(0.5, 0.5, 'Waveform data not found', 'Parent', ...
                ax1, 'HorizontalAlignment', 'center');
        end

        % Panel 2: ISI Histogram
        ax2 = subplot(6, 6, [3, 4]);
        if numel(spike_times) > 1
            isi = diff(spike_times) * 1000; % in ms
            histogram(ax2, isi, 'EdgeColor', 'k', 'FaceColor', ...
                [0.5 0.5 0.5]);
            set(ax2, 'XScale', 'log');
            title(ax2, 'ISI Histogram');
            xlabel(ax2, 'Inter-Spike Interval (ms, log scale)');
            ylabel(ax2, 'Count');
        else
            text(0.5, 0.5, 'Not enough spikes for ISI', 'Parent', ...
                ax2, 'HorizontalAlignment', 'center');
        end

        % Panel 3: PSTH for Tokens Task
        ax3 = subplot(6, 6, [5, 6]);
        if isfield(session_data, 'eventTimes') && isfield( ...
                session_data.eventTimes, 'pdsOutcomeOn')
            event_times = session_data.eventTimes.pdsOutcomeOn;
            win = [-0.5, 1.0];
            bin_size = 0.05;

            goodEvent = ~isnan(event_times) & event_times > 0;
            event_times_psth = session_data.eventTimes.trialBegin( ...
                goodEvent) + event_times(goodEvent);

            [~, psth, bin_centers] = alignAndBinSpikes(spike_times, ...
                event_times_psth, win(1), win(2), bin_size);
            n_trials = size(psth, 1);
            psth_rate = sum(psth, 1) / (n_trials * bin_size);

            if ~isempty(psth_rate)
                barStairsFill(bin_centers, psth_rate, zeros(size( ...
                    psth_rate)));
                title(ax3, 'PSTH around Outcome (Tokens)');
                xlabel(ax3, 'Time from Outcome (s)');
                ylabel(ax3, 'Firing Rate (Hz)');
                xline(ax3, 0, 'r--');
            else
                text(0.5, 0.5, 'No valid trials for PSTH', 'Parent', ...
                    ax3, 'HorizontalAlignment', 'center');
            end
        else
            text(0.5, 0.5, 'pdsOutcomeOn not found', 'Parent', ax3, ...
                'HorizontalAlignment', 'center');
        end

        % --- Middle Grid: Spatial Tuning ---
        gSac_trial_idx = session_data.trialInfo.taskCode == ...
            gSac_task_code;
        unique_thetas = unique(session_data.trialInfo.targetTheta( ...
            gSac_trial_idx));

        win = [-0.5, 1.0];
        bin_size = 0.05;

        saccade_axes = [];

        for i_theta = 1:numel(unique_thetas)
            current_theta = unique_thetas(i_theta);
            theta_trial_idx = gSac_trial_idx & ( ...
                session_data.trialInfo.targetTheta == current_theta);

            row_offset = (i_theta) * 6;

            % Left Column: Aligned to targetOn
            ax_target = subplot(6, 6, row_offset + [1:3]);
            event_times_target = session_data.eventTimes.targetOn( ...
                theta_trial_idx);
            [~, psth, bin_centers] = alignAndBinSpikes(spike_times, ...
                event_times_target, win(1), win(2), bin_size);
            psth_rate = sum(psth, 1) / (size(psth, 1) * bin_size);
            if ~isempty(psth_rate)
                barStairsFill(bin_centers, psth_rate, zeros(size( ...
                    psth_rate)));
                ylabel(ax_target, 'Firing Rate (Hz)');
                xline(ax_target, 0, 'r--');
            else
                text(0.5, 0.5, 'No data', 'Parent', ax_target, ...
                    'HorizontalAlignment', 'center');
            end
             if i_theta == 1
                title(ax_target, sprintf(['Aligned to Target ' ...
                    'On\nTheta: %d deg'], current_theta));
             else
                title(ax_target, sprintf('Theta: %d deg', ...
                    current_theta));
             end
             if i_theta < numel(unique_thetas)
                set(ax_target,'xticklabel',{[]})
             else
                xlabel(ax_target, 'Time from Target On (s)');
             end

            % Right Column: Aligned to saccadeOnset
            ax_saccade = subplot(6, 6, row_offset + [4:6]);
            event_times_saccade = session_data.eventTimes.saccadeOnset( ...
                theta_trial_idx);
            [~, psth, bin_centers] = alignAndBinSpikes(spike_times, ...
                event_times_saccade, win(1), win(2), bin_size);
            psth_rate = sum(psth, 1) / (size(psth, 1) * bin_size);
            if ~isempty(psth_rate)
                barStairsFill(bin_centers, psth_rate, zeros(size( ...
                    psth_rate)));
                xline(ax_saccade, 0, 'r--');
            else
                text(0.5, 0.5, 'No data', 'Parent', ax_saccade, ...
                    'HorizontalAlignment', 'center');
            end
            if i_theta == 1
                title(ax_saccade, 'Aligned to Saccade Onset');
            end
            if i_theta < numel(unique_thetas)
                set(ax_saccade,'xticklabel',{[]})
            else
                xlabel(ax_saccade, 'Time from Saccade Onset (s)');
            end
            saccade_axes = [saccade_axes, ax_saccade];
        end

        % Standardize Y-axis limits for saccade plots
        [~, ylims] = outerLims(saccade_axes);
        for i_ax = 1:numel(saccade_axes)
            set(saccade_axes(i_ax), 'YLim', ylims);
        end

        % --- Bottom Row: Summary Information ---
        ax_summary = subplot(6, 6, [31:36]);
        axis(ax_summary, 'off');

        screening_status = 'Not Selected';
        if numel(selected_neurons) >= i_cluster && selected_neurons( ...
                i_cluster)
            screening_status = 'Selected';
        end

        phy_quality = 'N/A';
        info_row = cluster_info.cluster_id == cluster_id;
        if any(info_row) && isfield(cluster_info, 'group')
            phy_quality = cluster_info.group(info_row);
        end

        baseline_fr = session_data.metrics.baseline_frs(i_cluster);
        wf_metrics = session_data.metrics.wf_metrics(i_cluster);
        wf_duration = sprintf('%.2f ms', wf_metrics.peak_trough_ms);

        summary_text = {
            sprintf('Cluster ID: %d', cluster_id), ...
            sprintf('Phy Quality: %s', phy_quality), ...
            sprintf('Screening Status: %s', screening_status), ...
            sprintf('Baseline FR: %.2f Hz', baseline_fr), ...
            sprintf('Waveform Duration: %s', wf_duration)
        };

        text(0.1, 0.5, summary_text, 'Parent', ax_summary, ...
            'VerticalAlignment', 'middle', 'FontSize', 12);
        title(ax_summary, 'Summary Information');
    else
        % --- Panel 1: Mean Waveform ---
        ax1 = subplot(2, 2, 1);
        if isfield(session_data.spikes, 'wfMeans') && numel( ...
                session_data.spikes.wfMeans) >= i_cluster
            plot(ax1, session_data.spikes.wfMeans{i_cluster}');
            title(ax1, 'Mean Waveform');
            xlabel(ax1, 'Samples');
            ylabel(ax1, 'Amplitude (uV)');
            axis(ax1, 'tight');
        else
            text(0.5, 0.5, 'Waveform data not found', 'Parent', ax1, ...
                'HorizontalAlignment', 'center');
        end

        % --- Panel 2: ISI Histogram ---
        ax2 = subplot(2, 2, 2);
        if numel(spike_times) > 1
            isi = diff(spike_times) * 1000; % in ms
            histogram(ax2, isi, 'EdgeColor', 'k', 'FaceColor', ...
                [0.5 0.5 0.5]);
            set(ax2, 'XScale', 'log');
            title(ax2, 'ISI Histogram');
            xlabel(ax2, 'Inter-Spike Interval (ms, log scale)');
            ylabel(ax2, 'Count');
        else
            text(0.5, 0.5, 'Not enough spikes for ISI', 'Parent', ...
                ax2, 'HorizontalAlignment', 'center');
        end

        % --- Panel 3: Basic PSTH ---
        ax3 = subplot(2, 2, 3);
        if isfield(session_data, 'eventTimes') && ...
                isfield(session_data.eventTimes, 'pdsOutcomeOn')
            event_times = session_data.eventTimes.pdsOutcomeOn;
            win = [-0.5, 1.0]; % Time window around event
            bin_size = 0.05; % 50 ms bins

            % get rid of NaN and -1 in 'pdsOutcomeOn'
            goodEvent = ~isnan(event_times) & event_times > 0;
            event_times_psth = session_data.eventTimes.trialBegin(...
                goodEvent) + event_times(goodEvent);

            [~, psth, bin_centers] = alignAndBinSpikes(spike_times, ...
                event_times_psth, win(1), win(2), bin_size);

            % Convert to firing rate
            n_trials = size(psth, 1);
            psth_rate = sum(psth, 1) / (n_trials * bin_size);

            if ~isempty(psth_rate)
                barStairsFill(bin_centers, psth_rate, zeros(size( ...
                    psth_rate)));
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
        ax4 = subplot(2, 2, 4);
        axis(ax4, 'off'); % Turn off axis for text

        screening_status = 'Not Selected';
        if numel(selected_neurons) >= i_cluster && ...
            selected_neurons(i_cluster)
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
    end

    % Appending only works with postscript files in MATLAB so we will
    % create one PDF per unit:
    pdfSave(output_filename, fig.Position(3:4)/72, fig);

    % give feedback:
    giveFeed(['Wrote PDF for unit ' sprintf('%03d',i_cluster)])
end

% Close the figure
close(fig);

fprintf('Diagnostic PDF saved to %s\n', output_filename);

end