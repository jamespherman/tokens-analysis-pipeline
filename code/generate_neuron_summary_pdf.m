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

% --- Setup and Data Extraction ---
codes = initCodes;
gSac_task_code = codes.uniqueTaskCode_gSac_4factors;
has_gSac_trials = any(session_data.trialInfo.taskCode == gSac_task_code);
has_gSac_events = isfield(session_data.eventTimes, 'targetOn') && ...
                  isfield(session_data.eventTimes, 'saccadeOnset');
do_spatial_tuning_plot = has_gSac_trials && has_gSac_events;

% Identify gSac_jph memory-guided saccade trials with valid reward times
gSac_jph_task_code = codes.uniqueTaskCode_gSac_jph;
has_gSac_jph_trials = any(session_data.trialInfo.taskCode == ...
    gSac_jph_task_code);
gSac_jph_memSac_trials = [];
if has_gSac_jph_trials
    if isfield(session_data.eventTimes, 'targetReillum') && ...
            isfield(session_data.eventTimes, 'pdsReward')
        gSac_jph_memSac_trials = find(session_data.trialInfo.taskCode ...
            == gSac_jph_task_code & ...
            ~isnan(session_data.eventTimes.targetReillum) & ...
            session_data.eventTimes.pdsReward > 0);
    else
        has_gSac_jph_trials = false;
    end
end

% We're going to lay the plots out in two columns, the 1st column is more
% diagnostic information than anything else, and the 2nd column is all
% PSTHs and rastergrams. We'll lay out the '1st column' assuming 6 rows of
% plots and 6 columns, using the 1st 3 columns of the top row for plotting
% and the reamining rows of the 1st the columns for text information. This
% is the '1st column' of diagnostic information I'm referring to above. The
% '2nd column' will be layed out assuming 10 rows of plots and 4 columns,
% and we'll use columns 3, and 4 for these plots. To do this we need to
% define the indexes for columns 3 and 4 under the 12 rows x 4 columns
% assumption so we can easily define the location of each of the PSTHs:
psthAndRasterPlotInd = reshape(1:(12*4), 4, 12)';
psthAndRasterPlotInd = psthAndRasterPlotInd(:, 3:4);

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
        % --- DYNAMICALLY DETERMINE SUBPLOT LAYOUT ---
        gSac_trial_idx = session_data.trialInfo.taskCode == gSac_task_code;
        unique_thetas = unique(session_data.trialInfo.targetTheta( ...
            gSac_trial_idx));
        n_thetas_gSac = numel(unique_thetas);

        unique_thetas_jph = [];
        n_thetas_jph = 0;
        if has_gSac_jph_trials && ~isempty(gSac_jph_memSac_trials)
            unique_thetas_jph = unique( ...
                session_data.trialInfo.targetTheta( ...
                gSac_jph_memSac_trials));
            n_thetas_jph = numel(unique_thetas_jph);
        end

        % --- Top Row: Original Diagnostics ---
        ax1 = subplot(6, 6, 1);
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

        ax2 = subplot(6, 6, 2);
        if numel(spike_times) > 1
            isi = diff(spike_times) * 1000; % in ms
            histogram(ax2, isi, 'EdgeColor', 'k', 'FaceColor', ...
                [0.5 0.5 0.5]);
            set(ax2, 'XScale', 'log');
            title(ax2, 'ISI Histogram');
            xlabel(ax2, 'Inter-Spike Interval (ms, log scale)');
            ylabel(ax2, 'Count');
        else
            text(0.5, 0.5, 'Not enough spikes for ISI', 'Parent', ax2, ...
                'HorizontalAlignment', 'center');
        end

        % This PSTH goes in the 2nd row, 1st column of the 'right column'
        % of panels, so  use psthAndRasterPlotInd(2,1) as the index.
        ax3 = subplot(12, 4, psthAndRasterPlotInd(2,1));
        if isfield(session_data, 'eventTimes') && isfield( ...
                session_data.eventTimes, 'pdsOutcomeOn')
            event_times = session_data.eventTimes.outcomeOn;
            win = [-0.5, 1.0];
            bin_size = 0.025;
            goodEvent = event_times > 0 & ...
                session_data.eventTimes.pdsReward > 0;
            event_times_psth = event_times(goodEvent);
            [~, psth, bin_centers] = alignAndBinSpikes(spike_times, ...
                event_times_psth, win(1), win(2), bin_size);
            n_trials = size(psth, 1);
            psth_rate = sum(psth, 1) / (n_trials * bin_size);
            if ~isempty(psth_rate)
                barStairsFill(bin_centers, psth_rate, zeros( ...
                    size(psth_rate)));
                
                xlabel(ax3, 'Time from Outcome (s)');
                ylabel(ax3, 'Firing Rate (Hz)');
                xline(ax3, 0, 'r--');
            else
                text(0.5, 0.5, 'No valid trials for PSTH', ...
                    'Parent', ax3, 'HorizontalAlignment', 'center');
            end
        else
            text(0.5, 0.5, 'pdsOutcomeOn not found', ...
                'Parent', ax3, 'HorizontalAlignment', 'center');
        end

        % the associated rastergram with 'ax3' goes directly above it:
        ax3r = subplot(12, 4, psthAndRasterPlotInd(1,1));
        imagesc(bin_centers, 1:n_trials, psth);
        colormap(flipud(bone(64)))
        set(ax3r, 'XLim', ax3.XLim, 'YLim', [0 n_trials+1])
        title(ax3r, 'PSTH around Outcome (Tokens)');

        % --- Middle Grid: Spatial Tuning (gSac_4factors) ---
        win = [-0.5, 1.0];
        bin_size = 0.025;
        for i_theta = 1:n_thetas_gSac
            current_theta = unique_thetas(i_theta);
            theta_trial_idx = gSac_trial_idx & ( ...
                session_data.trialInfo.targetTheta == current_theta);

            % the 1st theta target-onset-aligned data go in the 4th row,
            % 3rd column of plots. We abandon 'row_offset' and use our
            % index list instead, but with a sort of implicit row offset:
            ax_target = subplot(12, 4, ...
                psthAndRasterPlotInd(4 + (i_theta - 1) * 2, 1));
            event_times_target = session_data.eventTimes.targetOn(...
                theta_trial_idx);
            [~, psth, bin_centers] = alignAndBinSpikes(...
                spike_times, event_times_target, win(1), win(2), ...
                bin_size);
            psth_rate = sum(psth, 1) / (size(psth, 1) * bin_size);
            if ~isempty(psth_rate)
                barStairsFill(bin_centers, psth_rate, zeros( ...
                    size(psth_rate)));
                ylabel(ax_target, 'Firing Rate (Hz)');
                xline(ax_target, 0, 'r--');
            else
                text(0.5, 0.5, 'No data', 'Parent', ax_target, ...
                    'HorizontalAlignment', 'center');
            end

            % the rastergram associated with ax_target goes one row
            % above it:
            ax_target_r = subplot(12, 4, ...
                psthAndRasterPlotInd(4 + (i_theta - 1) * 2 - 1, 1));
            imagesc(bin_centers, 1:size(psth, 1), psth);
            colormap(flipud(bone(64)))
            set(ax_target_r, 'XLim', ax_target.XLim, 'YLim', ...
                [0 size(psth, 1)+1])
            if i_theta == 1
                title(ax_target_r, sprintf(['Aligned to Target ' ...
                    'On\nTheta: %d deg'], current_theta));
            else
                title(ax_target_r, sprintf('Theta: %d deg', current_theta));
            end
            if i_theta < n_thetas_gSac
                set(ax_target,'xticklabel',{[]})
            else
                xlabel(ax_target, 'Time from Target On (s)');
            end

            % the 1st theta target-onset-aligned data go in the 4th row,
            % 4th column of plots.
            ax_saccade = subplot(12, 4, ...
                psthAndRasterPlotInd(4 + (i_theta - 1) * 2, 2));
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
            

            % the rastergram associated with ax_saccade goes one row
            % above it:
            ax_saccade_r = subplot(12, 4, ...
                psthAndRasterPlotInd(4 + (i_theta - 1) * 2 - 1, 2));
            imagesc(bin_centers, 1:size(psth, 1), psth);
            colormap(flipud(bone(64)))
            set(ax_saccade_r, 'XLim', ax_target.XLim, 'YLim', ...
                [0 size(psth, 1)+1])
            if i_theta == 1
                title(ax_saccade, 'Aligned to Saccade Onset');
            end
            if i_theta < n_thetas_gSac
                set(ax_saccade,'xticklabel',{[]})
            else
                xlabel(ax_saccade, 'Time from Saccade Onset (s)');
            end
        end

        % --- Middle Grid: Spatial Tuning (gSac_jph) ---
        if n_thetas_jph > 0
            for i_theta_jph = 1:n_thetas_jph
                current_theta = unique_thetas_jph(i_theta_jph);
                theta_trial_idx_in_mem_sac = ...
                    session_data.trialInfo.targetTheta(...
                    gSac_jph_memSac_trials) == current_theta;
                final_trial_indices = gSac_jph_memSac_trials(...
                    theta_trial_idx_in_mem_sac);

                % gSac_jph target PSTH axes go at the bottom, 3rd column
                ax_target_jph = subplot(12, 4, ...
                    psthAndRasterPlotInd(12, 1));
                event_times_target = session_data.eventTimes.targetOn(...
                    final_trial_indices);
                [~, psth, bin_centers] = alignAndBinSpikes(spike_times, ...
                    event_times_target, win(1), win(2), bin_size);
                psth_rate = sum(psth, 1) / (size(psth, 1) * bin_size);
                if ~isempty(psth_rate)
                    barStairsFill(bin_centers, psth_rate, zeros(...
                        size(psth_rate)));
                    ylabel(ax_target_jph, 'Firing Rate (Hz)');
                    xline(ax_target_jph, 0, 'r--');
                else
                    text(0.5, 0.5, 'No data', 'Parent', ax_target_jph, ...
                        'HorizontalAlignment', 'center');
                end
                
                if i_theta_jph < n_thetas_jph
                    set(ax_target_jph,'xticklabel',{[]})
                else
                    xlabel(ax_target_jph, 'Time from Target On (s)');
                end

                % gSac_jph target rastergram axes go above the PSTH axes
                ax_target_jph_r = subplot(12, 4, ...
                    psthAndRasterPlotInd(11, 1));
                imagesc(bin_centers, 1:size(psth, 1), psth);
                colormap(flipud(bone(64)))
                set(ax_target_jph_r, 'XLim', ax_target_jph.XLim, 'YLim', ...
                    [0 size(psth, 1)+1])
                title(ax_target_jph, {sprintf(['PSTH around Target On ' ...
                    '(gSac_jph Mem)']), sprintf('Theta: %d deg', ...
                    current_theta)});

                % gSac_jph saccade PSTH axes go in the bottom row (12) 4th
                % column:
                ax_saccade_jph = subplot(12, 4, ...
                    psthAndRasterPlotInd(12, 2));
                event_times_saccade = ...
                    session_data.eventTimes.saccadeOnset(...
                    final_trial_indices);
                [~, psth, bin_centers] = alignAndBinSpikes(...
                    spike_times, event_times_saccade, win(1), win(2), ...
                    bin_size);
                psth_rate = sum(psth, 1) / (size(psth, 1) * bin_size);
                if ~isempty(psth_rate)
                    barStairsFill(bin_centers, psth_rate, zeros(...
                        size(psth_rate)));
                    xline(ax_saccade_jph, 0, 'r--');
                else
                    text(0.5, 0.5, 'No data', 'Parent', ax_saccade_jph, ...
                        'HorizontalAlignment', 'center');
                end
                
                if i_theta_jph < n_thetas_jph
                    set(ax_saccade_jph,'xticklabel',{[]})
                else
                    xlabel(ax_saccade_jph, 'Time from Saccade Onset (s)');
                end

                % gSac_jph saccade rastergram axes go above the PSTH axes
                ax_saccade_jph_r = subplot(12, 4, ...
                    psthAndRasterPlotInd(11, 2));
                imagesc(bin_centers, 1:size(psth, 1), psth);
                colormap(flipud(bone(64)))
                set(ax_saccade_jph_r, 'XLim', ax_saccade_jph.XLim, 'YLim', ...
                    [0 size(psth, 1)+1])
                if i_theta_jph == 1
                    title(ax_saccade_jph, ...
                        'Aligned to Saccade Onset (gSac_jph Mem)');
                end
            end
        end

        % --- Bottom Row: Summary Information ---
        ax_summary = subplot(6, 6, [25:27 31:33]);
        axis(ax_summary, 'off');

        screening_status = 'Not Selected';
        if numel(selected_neurons) >= i_cluster && ...
                selected_neurons(i_cluster)
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
        ax1 = subplot(6, 6, 1);
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
        ax2 = subplot(6, 6, 2);
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
        ax3 = subplot(12, 4, [11 15]);
        if isfield(session_data, 'eventTimes') && ...
                isfield(session_data.eventTimes, 'pdsOutcomeOn')
            event_times = session_data.eventTimes.pdsOutcomeOn;
            win = [-0.5, 1.0]; % Time window around event
            bin_size = 0.025; % 50 ms bins

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

        % the associated rastergram with 'ax3' goes directly above it:
        ax3r = subplot(12, 4, [3 7]);
        imagesc(bin_centers, 1:n_trials, psth);
        colormap(flipud(bone(64)))
        set(ax3r, 'XLim', ax3.XLim, 'YLim', [0 n_trials+1])

        % --- Panel 4: Summary Information ---
        ax4 = subplot(6, 6, [5 11 17 6 12 18]);
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