function generate_neuron_summary_pdf(session_data, ...
    selected_neurons, unique_id)
% GENERATE_NEURON_SUMMARY_PDF Generates a multi-page PDF with
% diagnostic plots for each neuron.
%
%   GENERATE_NEURON_SUMMARY_PDF(session_data, selected_neurons, ...
%   unique_id)
%   generates a multi-page PDF file where each page contains a summary
%   of a single neuron, conforming to the standardized `session_data`
%   structure.
%
%   Inputs:
%       session_data     - A struct for a single session, containing
%                          fields like 'spikes' and 'cluster_info'.
%       selected_neurons - A logical vector (nClusters x 1) indicating
%                          which neurons were selected by a screening
%                          process.
%       unique_id        - A string used to name the output PDF file.
%
%   Output:
%       The function saves a PDF file named
%       '[unique_id]_neuron_diagnostics.pdf' in the 'figures/'
%       directory.

%% Setup Paths
% Add the 'utils' directory to the path so that helper functions can
% be found.
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Constants
% Define constants for plot layout and analysis parameters.
N_ROWS = 8;
N_COLS = 6;
PSTH_WINDOW = [-0.5, 1.0]; % 1.5s window for PSTHs
PSTH_BIN_SIZE = 0.025;     % 25ms bin size

%% Configuration
% In-line function to report timing
giveFeed = @(x)disp([num2str(toc) 's - ' x]);

% Define output directory and filename
project_root = fullfile(findOneDrive, 'Code', ...
    'tokens-analysis-pipeline');
output_dir = fullfile(project_root, 'figures');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% --- Setup and Data Extraction ---
cluster_info = session_data.spikes.cluster_info;
cluster_ids = cluster_info.cluster_id;
nClusters = height(cluster_info.cluster_id);
all_spike_times = session_data.spikes.times;
all_spike_clusters = session_data.spikes.clusters;

% --- Check for required data ---
codes = initCodes;
gSac_task_code = codes.uniqueTaskCode_gSac_4factors;
has_gSac_trials = any(session_data.trialInfo.taskCode == ...
    gSac_task_code);
has_gSac_events = isfield(session_data.eventTimes, 'targetOn') && ...
                  isfield(session_data.eventTimes, 'saccadeOnset');
do_spatial_tuning_plot = has_gSac_trials && has_gSac_events;

% Identify gSac_jph memory-guided saccade trials with valid reward
% times
gSac_jph_task_code = codes.uniqueTaskCode_gSac_jph;
has_gSac_jph_trials = any(session_data.trialInfo.taskCode == ...
    gSac_jph_task_code);
gSac_jph_memSac_trials = [];
if has_gSac_jph_trials
    has_required_fields = isfield(session_data.eventTimes, ...
        'targetReillum') && isfield(session_data.eventTimes, ...
        'pdsReward');
    if has_required_fields
        gSac_jph_memSac_trials = find(...
            session_data.trialInfo.taskCode == gSac_jph_task_code & ...
            ~isnan(session_data.eventTimes.targetReillum) & ...
            session_data.eventTimes.pdsReward > 0);
    else
        has_gSac_jph_trials = false;
    end
end

% Loop through each cluster to create a page of plots
for i_cluster = 1:nClusters

    % Create a new figure for the PDF
    fig = figure('Color', 'w', 'MenuBar', 'None', 'ToolBar', 'None');

    % define CLUSTER ID:
    cluster_id = cluster_ids(i_cluster);
    
    % define output file name:
    output_filename = fullfile(output_dir, ...
        [unique_id '_' sprintf('%03d',i_cluster) ...
        '_neuron_diagnostics.pdf']);

    % Add a title for the page
    sgtitle(sprintf('Neuron Diagnostic Summary: Cluster %d', cluster_id));

    spike_times = all_spike_times(all_spike_clusters == cluster_id);

    if do_spatial_tuning_plot
        %% --- Spatial Tuning Plots ---

        % --- Data Setup for Spatial Tuning ---
        gSac_trial_idx = session_data.trialInfo.taskCode == ...
            gSac_task_code;
        unique_thetas = unique(session_data.trialInfo.targetTheta( ...
            gSac_trial_idx));
        n_thetas_gSac = numel(unique_thetas);

        unique_thetas_jph = [];
        if has_gSac_jph_trials && ~isempty(gSac_jph_memSac_trials)
            unique_thetas_jph = unique( ...
                session_data.trialInfo.targetTheta( ...
                gSac_jph_memSac_trials));
        end
        n_thetas_jph = numel(unique_thetas_jph);

        % --- Column 1-2: Basic Diagnostics ---
        % Waveform Plot
        ax_wf = mySubPlot([N_ROWS, N_COLS, 1]);
        if isfield(session_data.spikes, 'wfMeans') && ...
                numel(session_data.spikes.wfMeans) >= i_cluster
            plot(ax_wf, session_data.spikes.wfMeans{i_cluster}');
            title(ax_wf, 'Mean Waveform');
            xlabel(ax_wf, 'Samples');
            ylabel(ax_wf, 'Amplitude (uV)');
            axis(ax_wf, 'tight');
        else
            text(0.5, 0.5, 'No waveform', 'Parent', ax_wf, ...
                'HorizontalAlignment', 'center');
        end

        % ISI Histogram
        ax_isi = mySubPlot([N_ROWS, N_COLS, 2]);
        if numel(spike_times) > 1
            isi = diff(spike_times) * 1000; % ms
            histogram(ax_isi, isi, 'EdgeColor', 'k', ...
                'FaceColor', [0.5 0.5 0.5]);
            set(ax_isi, 'XScale', 'log');
            title(ax_isi, 'ISI Histogram');
            xlabel(ax_isi, 'ISI (ms, log)');
            ylabel(ax_isi, 'Count');
        else
            text(0.5, 0.5, 'No ISI', 'Parent', ax_isi, ...
                'HorizontalAlignment', 'center');
        end

        % --- Column 3-4: PSTH for Tokens Task Outcome ---
        ax_raster = mySubPlot([N_ROWS, N_COLS, 3]);
        ax_psth = mySubPlot([N_ROWS, N_COLS, 3 + N_COLS]);
        if isfield(session_data.eventTimes, 'outcomeOn')
            event_times = session_data.eventTimes.outcomeOn;
            valid_trials = event_times > 0 & ...
                session_data.eventTimes.pdsReward > 0;
            plot_psth_and_raster(ax_raster, ax_psth, spike_times, ...
                event_times(valid_trials), PSTH_WINDOW, ...
                PSTH_BIN_SIZE, 'Outcome (Tokens)', ...
                'Time from Outcome (s)', true);
        else
            text(0.5, 0.5, 'No outcome event', 'Parent', ax_raster, ...
                'HorizontalAlignment', 'center');
        end

        % --- Column 5-6: PSTHs for gSac Tasks ---
        % gSac_4factors task
        for i_theta = 1:n_thetas_gSac
            row_offset = (i_theta - 1) * 2;
            current_theta = unique_thetas(i_theta);
            theta_trials = gSac_trial_idx & ...
                (session_data.trialInfo.targetTheta == current_theta);

            % Aligned to Target On
            ax_r = mySubPlot([N_ROWS, N_COLS, 5 + row_offset]);
            ax_p = mySubPlot([N_ROWS, N_COLS, 5 + N_COLS + row_offset]);
            et = session_data.eventTimes.targetOn(theta_trials);
            plot_psth_and_raster(ax_r, ax_p, spike_times, et, ...
                PSTH_WINDOW, PSTH_BIN_SIZE, ...
                sprintf('Target On (gSac): %d deg', current_theta), ...
                'Time from Target On (s)', i_theta == 1);
            if i_theta < n_thetas_gSac, set(ax_p, 'xticklabel', {[]}); end

            % Aligned to Saccade Onset
            ax_r = mySubPlot([N_ROWS, N_COLS, 6 + row_offset]);
            ax_p = mySubPlot([N_ROWS, N_COLS, 6 + N_COLS + row_offset]);
            et = session_data.eventTimes.saccadeOnset(theta_trials);
            plot_psth_and_raster(ax_r, ax_p, spike_times, et, ...
                PSTH_WINDOW, PSTH_BIN_SIZE, 'Saccade Onset (gSac)',...
                'Time from Saccade Onset (s)', false);
            if i_theta < n_thetas_gSac, set(ax_p, 'xticklabel', {[]}); end
        end

        % gSac_jph task
        if n_thetas_jph > 0
            row_offset = n_thetas_gSac * 2;
            i_theta_jph = 1; % For now, only plot the first one
            current_theta = unique_thetas_jph(i_theta_jph);
            theta_trials_in_mem_sac = ...
                session_data.trialInfo.targetTheta( ...
                gSac_jph_memSac_trials) == current_theta;
            final_indices = gSac_jph_memSac_trials(...
                theta_trials_in_mem_sac);

            % Aligned to Target On
            ax_r = mySubPlot([N_ROWS, N_COLS, 5 + row_offset]);
            ax_p = mySubPlot([N_ROWS, N_COLS, 5 + N_COLS + row_offset]);
            et = session_data.eventTimes.targetOn(final_indices);
            plot_psth_and_raster(ax_r, ax_p, spike_times, et, ...
                PSTH_WINDOW, PSTH_BIN_SIZE, ...
                sprintf('Target On (jph): %d deg', current_theta), ...
                'Time from Target On (s)', true);

            % Aligned to Saccade Onset
            ax_r = mySubPlot([N_ROWS, N_COLS, 6 + row_offset]);
            ax_p = mySubPlot([N_ROWS, N_COLS, 6 + N_COLS + row_offset]);
            et = session_data.eventTimes.saccadeOnset(final_indices);
            plot_psth_and_raster(ax_r, ax_p, spike_times, et, ...
                PSTH_WINDOW, PSTH_BIN_SIZE, ...
                'Saccade Onset (jph)', ...
                'Time from Saccade Onset (s)', false);
        end
    else
        %% --- Basic Plots (No Spatial Tuning) ---
        % Waveform Plot
        ax_wf = mySubPlot([N_ROWS, N_COLS, 1]);
        if isfield(session_data.spikes, 'wfMeans') && ...
                numel(session_data.spikes.wfMeans) >= i_cluster
            plot(ax_wf, session_data.spikes.wfMeans{i_cluster}');
            title(ax_wf, 'Mean Waveform');
            xlabel(ax_wf, 'Samples');
            ylabel(ax_wf, 'Amplitude (uV)');
            axis(ax_wf, 'tight');
        else
            text(0.5, 0.5, 'No waveform', 'Parent', ax_wf, ...
                'HorizontalAlignment', 'center');
        end

        % ISI Histogram
        ax_isi = mySubPlot([N_ROWS, N_COLS, 2]);
        if numel(spike_times) > 1
            isi = diff(spike_times) * 1000; % ms
            histogram(ax_isi, isi, 'EdgeColor', 'k', 'FaceColor', ...
                [0.5 0.5 0.5]);
            set(ax_isi, 'XScale', 'log');
            title(ax_isi, 'ISI Histogram');
            xlabel(ax_isi, 'ISI (ms, log)');
            ylabel(ax_isi, 'Count');
        else
            text(0.5, 0.5, 'No ISI', 'Parent', ax_isi, ...
                'HorizontalAlignment', 'center');
        end

        % PSTH for Tokens Task Outcome
        ax_raster = mySubPlot([N_ROWS, N_COLS, 3]);
        ax_psth = mySubPlot([N_ROWS, N_COLS, 3 + N_COLS]);
        if isfield(session_data.eventTimes, 'outcomeOn')
            event_times = session_data.eventTimes.outcomeOn;
            valid_trials = event_times > 0 & ...
                session_data.eventTimes.pdsReward > 0;
            plot_psth_and_raster(ax_raster, ax_psth, spike_times, ...
                event_times(valid_trials), PSTH_WINDOW, ...
                PSTH_BIN_SIZE, 'Outcome (Tokens)', ...
                'Time from Outcome (s)', true);
        else
            text(0.5, 0.5, 'No outcome event', 'Parent', ax_raster, ...
                'HorizontalAlignment', 'center');
        end
    end

    %% --- Summary Information ---
    ax_summary = mySubPlot([N_ROWS, N_COLS, N_COLS * (N_ROWS - 1) + 1]);
    axis(ax_summary, 'off');

    screening_status = 'Not Selected';
    if numel(selected_neurons) >= i_cluster && selected_neurons(i_cluster)
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
        'VerticalAlignment', 'middle', 'FontSize', 10);
    title(ax_summary, 'Summary Information');

    % modify appearance:
    set(findall(fig, 'Type', 'Axes'), 'XColor', 'k', 'YColor', 'k', ...
        'TickDir', 'Out', 'LineWidth', 1, 'Color', 'w', 'Box', 'Off');
    set(findall(fig, 'Type', 'Text'), 'Color', 0.2*[1 1 1], ...
        'FontSize', 8);
    set(fig, 'Position', [150 50 750 1000]);
    drawnow;

    % Appending only works with postscript files in MATLAB so we will
    % create one PDF per unit:
    pdfSave(output_filename, fig.Position(3:4)/72, fig);

    % give feedback:
    giveFeed(['Wrote PDF for unit ' sprintf('%03d',i_cluster)])

    % Close the figure
    close(fig);
end

fprintf('Diagnostic PDF saved to %s\n', output_filename);

end

%% Local Helper Functions

function plot_psth_and_raster(ax_raster, ax_psth, spike_times, ...
    event_times, psth_win, bin_size, title_text, xlabel_text, ...
    show_y_label)
% PLOT_PSTH_AND_RASTER Plots a raster and its corresponding PSTH.
%
%   Inputs:
%       ax_raster      - Axes handle for the raster plot.
%       ax_psth        - Axes handle for the PSTH plot.
%       spike_times    - Vector of spike timestamps.
%       event_times    - Vector of event timestamps to align to.
%       psth_win       - 2-element vector specifying the time window
%                        [start, end] around the event.
%       bin_size       - Width of the bins for the PSTH in seconds.
%       title_text     - String for the title of the raster plot.
%       xlabel_text    - String for the xlabel of the PSTH plot.
%       show_y_label   - Boolean to control visibility of the PSTH
%                        y-label.

% Calculate PSTH
[~, psth, bin_centers] = alignAndBinSpikes(spike_times, ...
    event_times, psth_win(1), psth_win(2), bin_size);
n_trials = size(psth, 1);
psth_rate = sum(psth, 1) / (n_trials * bin_size);

% Plot Raster
imagesc(ax_raster, bin_centers, 1:n_trials, psth);
colormap(ax_raster, flipud(bone(64)));
set(ax_raster, 'YDir', 'normal');
title(ax_raster, title_text);
set(ax_raster, 'xticklabel', {[]}); % No x-labels on raster

% Plot PSTH
if ~isempty(psth_rate)
    axes(ax_psth); % Select the correct axes
    barStairsFill(bin_centers, psth_rate, zeros(size(psth_rate)));
    xline(ax_psth, 0, 'r--');
    xlabel(ax_psth, xlabel_text);
    if show_y_label
        ylabel(ax_psth, 'Firing Rate (Hz)');
    else
        set(ax_psth, 'yticklabel', {[]});
    end
else
    text(0.5, 0.5, 'No data', 'Parent', ax_psth, ...
        'HorizontalAlignment', 'center');
end

% Link axes and set limits
linkaxes([ax_raster, ax_psth], 'x');
set(ax_raster, 'XLim', psth_win, 'YLim', [0 n_trials + 1]);

end