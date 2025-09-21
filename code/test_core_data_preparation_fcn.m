function test_core_data_preparation_fcn(unique_id)
%% test_core_data_preparation_fcn.m
%
%   This function tests the refactored `prepare_core_data` module. It loads
%   a sample session, runs the data preparation pipeline, and generates a
%   figure with diagnostic plots to verify the output for all three key
%   alignment events: CUE_ON, outcomeOn, and reward.
%
%   Args:
%       unique_id (char): The unique identifier for the session to be tested.
%
% Author: Jules
% Date: 2025-09-21

%% Setup
% Start timer
tic;

% In-line function to report timing
giveFeed = @(x)disp([num2str(toc) 's - ' x]);

% Add utility functions to the MATLAB path
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Intro
giveFeed(sprintf('Testing diagnostic workflow for session: %s', ...
    unique_id));

%% Load Manifest
giveFeed('Step 1: Loading session manifest...');
project_root = fullfile(findOneDrive, 'Code', 'tokens-analysis-pipeline');
manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');
manifest = readtable(manifest_path);

% Find the row for the current session
session_idx = find(strcmp(manifest.unique_id, unique_id));
if isempty(session_idx)
    error('test_core_data_preparation:sessionNotFound', ...
          'Session with unique_id "%s" not found in the manifest.', ...
          unique_id);
end

%% Load Data
% Define the project root and construct the path to a sample data file.
giveFeed('Step 2: Loading session data...');
one_drive_path = findOneDrive;
data_file = fullfile(one_drive_path, ...
    'Neuronal Data Analysis', unique_id, ...
    [unique_id '_session_data.mat']);

% Check if the sample file exists
if ~exist(data_file, 'file')
    error('test_core_data_preparation:sampleNotFound', ...
          ['Sample data file not found. Please specify a valid path ' ...
          'to a \n session_data.mat file in the "Load Data" section of ' ...
          'this script.']);
end

fprintf('Loading session data from: %s\n', data_file);
load(data_file, 'session_data');

%% Step 2: Calculate and Store Diagnostic Metrics (if needed)
if ~isfield(session_data, 'metrics')
giveFeed('Step 2: Calculating and storing diagnostic metrics...');

% Calculate baseline firing rates and add to session_data
session_data.metrics.baseline_frs = calculate_baseline_fr(session_data);

% Calculate waveform metrics and add to session_data
nClusters = height(session_data.spikes.cluster_info.cluster_id);
for i_cluster = 1:nClusters

    % get current unit's multi-channel waveform:
    mean_waveform = session_data.spikes.wfMeans{i_cluster};

    % find channel with max variance:
    [~,max_var_chan] = max(var(mean_waveform,[],2));

    % Note: Assuming a sampling rate of 30000 Hz
    session_data.metrics.wf_metrics(i_cluster, 1) = ...
        calculate_waveform_metrics(mean_waveform(max_var_chan,:), 30000);
end
giveFeed('Diagnostic metrics calculated and stored.');
end

%% Run neuron screening
giveFeed('Step 3: Running neuron screening...');
if strcmp(manifest.screening_status{session_idx}, 'complete')
    giveFeed(['Screening status is ''complete''. ' ...
        'Loading results from session_data...']);

    % Load pre-computed results from the session_data struct
    selected_neurons = session_data.analysis.selected_neurons;
    if isfield(session_data.analysis, 'scSide')
        scSide = session_data.analysis.scSide;
    end

else
    giveFeed('Screening status is ''pending''. Running screening functions...');
    selected_neurons = []; % Initialize empty
    session_data.metadata.unique_id = unique_id;
    if contains(unique_id, 'SNc')
        giveFeed('Session is SNc type. Running screen_da_neurons...');
        selected_neurons = screen_da_neurons(session_data, unique_id);
        giveFeed('DA neuron screening complete.');

    elseif contains(unique_id, 'SC')
        giveFeed('Session is SC type. Running screen_sc_neurons...');

        % screen_sc_neurons will now determine scSide and calculate significance
        [selected_neurons, sig_epoch_comp, scSide] = screen_sc_neurons(session_data);

        % Store results in a new 'analysis' field
        session_data.analysis.scSide = scSide;
        session_data.analysis.sig_epoch_comparison = sig_epoch_comp;

        giveFeed(sprintf('SC neuron screening complete. Determined scSide: %s.', scSide));
    else
        error('Unknown session type. Cannot determine which screening function to run.');
    end
    
    % Add selected_neurons to the analysis field
    session_data.analysis.selected_neurons = selected_neurons;

    % Save the updated session_data struct
    giveFeed('Saving screening results to session_data.mat...');
    save(data_file, 'session_data', '-v7.3');

    % Update manifest status in memory
    manifest.screening_status{session_idx} = 'complete';
    giveFeed('Manifest screening_status updated to ''complete''.');
end


%% Run Core Data Preparation
% This is the main function being tested
giveFeed('Step 4: Preparing core data...');
if strcmp(manifest.dataprep_status{session_idx}, 'complete')
    giveFeed('Data prep status is ''complete''. Loading core_data from session_data...');
    core_data = session_data.analysis.core_data;
else
    giveFeed('Data prep status is ''pending''. Running prepare_core_data...');
    core_data = prepare_core_data(session_data, selected_neurons);

    % Save the core_data struct into the session_data.analysis field
    session_data.analysis.core_data = core_data;

    % Save the updated session_data struct
    giveFeed('Saving core_data to session_data.mat...');
    save(data_file, 'session_data', '-v7.3');

    % Update manifest status in memory
    manifest.dataprep_status{session_idx} = 'complete';
    giveFeed('Manifest dataprep_status updated to ''complete''.');
end

%% Define Task Conditions
% Create logical masks for different trial conditions
giveFeed('Step 5: Defining task conditions...');
[conditions, is_av_session] = define_task_conditions(session_data.trialInfo, ...
    session_data.eventTimes, session_data.metadata.unique_id);

%% PLOTTING CONFIGURATION
giveFeed('Step 6: Configuring plots...');

% --- Define Colors ---
colors = struct();
colors.normal_dist = [0.2, 0.2, 0.8]; % Blue
colors.uniform_dist = [0.8, 0.2, 0.2]; % Red
colors.norm_common = [0.5, 0.5, 0.5]; % Gray
colors.norm_rare_high = [0.9, 0.6, 0];   % Orange
colors.flicker_surprising = [0.8, 0.2, 0.8]; % Magenta
colors.flicker_certain = [0.2, 0.8, 0.2]; % Green

% --- Pre-calculate data for all plots ---
if any(selected_neurons)
    % --- RPE Data ---
    % CUE_ON
    rates = core_data.spikes.CUE_ON.rates;
    [grand_mean_norm, psth_norm] = calculate_mean_psth(rates, conditions.is_normal_dist);
    [grand_mean_unif, psth_unif] = calculate_mean_psth(rates, conditions.is_uniform_dist);
    % outcomeOn
    rates = core_data.spikes.outcomeOn.rates;
    [grand_mean_common, psth_common] = calculate_mean_psth(rates, conditions.is_norm_common);
    [grand_mean_rare, psth_rare] = calculate_mean_psth(rates, conditions.is_norm_rare_high);
    % reward
    rates = core_data.spikes.reward.rates;
    [grand_mean_common_rew, psth_common_rew] = calculate_mean_psth(rates, conditions.is_norm_common);
    [grand_mean_rare_rew, psth_rare_rew] = calculate_mean_psth(rates, conditions.is_norm_rare_high);

    % --- SPE Data (if applicable) ---
    if is_av_session
        % CUE_ON
        rates = core_data.spikes.CUE_ON.rates;
        [grand_mean_surprising_cue, psth_surprising_cue] = calculate_mean_psth(rates, conditions.is_flicker_surprising);
        [grand_mean_certain_cue, psth_certain_cue] = calculate_mean_psth(rates, conditions.is_flicker_certain);
        % outcomeOn
        rates = core_data.spikes.outcomeOn.rates;
        [grand_mean_surprising_outcome, psth_surprising_outcome] = calculate_mean_psth(rates, conditions.is_flicker_surprising);
        [grand_mean_certain_outcome, psth_certain_outcome] = calculate_mean_psth(rates, conditions.is_flicker_certain);
        % reward
        rates = core_data.spikes.reward.rates;
        [grand_mean_surprising_reward, psth_surprising_reward] = calculate_mean_psth(rates, conditions.is_flicker_surprising);
        [grand_mean_certain_reward, psth_certain_reward] = calculate_mean_psth(rates, conditions.is_flicker_certain);
    end

    % --- Determine common color scale for heatmaps ---
    all_psth = {psth_norm, psth_unif, psth_common, psth_rare, psth_common_rew, psth_rare_rew};
    if is_av_session
        all_psth = [all_psth, {psth_surprising_cue, psth_certain_cue, ...
            psth_surprising_outcome, psth_certain_outcome, ...
            psth_surprising_reward, psth_certain_reward}];
    end
    all_maxes = cellfun(@(x) max(x, [], 'all'), all_psth);
    max_rate = max(all_maxes);
    c_lims = [0, max_rate];
end

% --- Master Figure Configuration ---
figures = struct();

% --- RPE Figure Definition ---
figures.RPE.title = sprintf('Core Data Verification (RPE): %s', unique_id);
figures.RPE.fileName = fullfile(project_root, 'figures', [unique_id '_RPE_verification.pdf']);
figures.RPE.panel_config = [];

if any(selected_neurons)
    align_events = {'CUE_ON', 'outcomeOn', 'reward'};
    psth_data_rpe = {{psth_norm, psth_unif}, {psth_common, psth_rare}, {psth_common_rew, psth_rare_rew}};
    grand_mean_data_rpe = {{grand_mean_norm, grand_mean_unif}, {grand_mean_common, grand_mean_rare}, {grand_mean_common_rew, grand_mean_rare_rew}};
    pupil_conditions_rpe = {{conditions.is_normal_dist, conditions.is_uniform_dist}, {conditions.is_norm_common, conditions.is_norm_rare_high}, {conditions.is_norm_common, conditions.is_norm_rare_high}};
    color_sets_rpe = {{colors.normal_dist, colors.uniform_dist}, {colors.norm_common, colors.norm_rare_high}, {colors.norm_common, colors.norm_rare_high}};
    legend_sets_rpe = {{'Normal Dist', 'Uniform Dist'}, {'Common', 'Rare High'}, {'Common', 'Rare High'}};

    p = 1;
    % Row 1: Heatmaps (cond 1)
    for i_align = 1:numel(align_events)
        align_event = align_events{i_align};
        figures.RPE.panel_config(p).dataType = 'Heatmap';
        figures.RPE.panel_config(p).alignEvent = align_event;
        figures.RPE.panel_config(p).xdata = {core_data.spikes.(align_event).time_vector};
        figures.RPE.panel_config(p).ydata = {psth_data_rpe{i_align}{1}};
        figures.RPE.panel_config(p).c_lims = c_lims;
        p = p+1;
    end
    % Row 2: Heatmaps (cond 2)
    for i_align = 1:numel(align_events)
        align_event = align_events{i_align};
        figures.RPE.panel_config(p).dataType = 'Heatmap';
        figures.RPE.panel_config(p).alignEvent = align_event;
        figures.RPE.panel_config(p).xdata = {core_data.spikes.(align_event).time_vector};
        figures.RPE.panel_config(p).ydata = {psth_data_rpe{i_align}{2}};
        figures.RPE.panel_config(p).c_lims = c_lims;
        p = p+1;
    end
    % Row 3: PSTHs
    for i_align = 1:numel(align_events)
        align_event = align_events{i_align};
        figures.RPE.panel_config(p).dataType = 'PSTH';
        figures.RPE.panel_config(p).alignEvent = align_event;
        figures.RPE.panel_config(p).xdata = {core_data.spikes.(align_event).time_vector, core_data.spikes.(align_event).time_vector};
        figures.RPE.panel_config(p).ydata = grand_mean_data_rpe{i_align};
        figures.RPE.panel_config(p).colors = cell2mat(color_sets_rpe{i_align}');
        figures.RPE.panel_config(p).legendStrings = legend_sets_rpe{i_align};
        p = p+1;
    end
    % Row 4: Pupil
    for i_align = 1:numel(align_events)
        align_event = align_events{i_align};
        figures.RPE.panel_config(p).dataType = 'Pupil';
        figures.RPE.panel_config(p).alignEvent = align_event;
        figures.RPE.panel_config(p).xdata = {core_data.pupil.(align_event).time_vector, core_data.pupil.(align_event).time_vector};
        mean_trace1 = mean(core_data.pupil.(align_event).traces(pupil_conditions_rpe{i_align}{1}, :), 1, 'omitnan');
        mean_trace2 = mean(core_data.pupil.(align_event).traces(pupil_conditions_rpe{i_align}{2}, :), 1, 'omitnan');
        figures.RPE.panel_config(p).ydata = {mean_trace1, mean_trace2};
        figures.RPE.panel_config(p).colors = cell2mat(color_sets_rpe{i_align}');
        p = p+1;
    end
end

% --- SPE Figure Definition ---
if is_av_session
    figures.SPE.title = sprintf('Core Data Verification (SPE): %s', unique_id);
    figures.SPE.fileName = fullfile(project_root, 'figures', [unique_id '_SPE_verification.pdf']);
    figures.SPE.panel_config = [];

    if any(selected_neurons)
        align_events = {'CUE_ON', 'outcomeOn', 'reward'};
        psth_data_spe = {{psth_surprising_cue, psth_certain_cue}, {psth_surprising_outcome, psth_certain_outcome}, {psth_surprising_reward, psth_certain_reward}};
        grand_mean_data_spe = {{grand_mean_surprising_cue, grand_mean_certain_cue}, {grand_mean_surprising_outcome, grand_mean_certain_outcome}, {grand_mean_surprising_reward, grand_mean_certain_reward}};
        pupil_conditions_spe = {{conditions.is_flicker_surprising, conditions.is_flicker_certain}, {conditions.is_flicker_surprising, conditions.is_flicker_certain}, {conditions.is_flicker_surprising, conditions.is_flicker_certain}};
        color_sets_spe = {{colors.flicker_surprising, colors.flicker_certain}, {colors.flicker_surprising, colors.flicker_certain}, {colors.flicker_surprising, colors.flicker_certain}};
        legend_sets_spe = {{'Surprising', 'Certain'}, {'Surprising', 'Certain'}, {'Surprising', 'Certain'}};

        p = 1;
        % Row 1: Heatmaps (cond 1)
        for i_align = 1:numel(align_events)
            align_event = align_events{i_align};
            figures.SPE.panel_config(p).dataType = 'Heatmap';
            figures.SPE.panel_config(p).alignEvent = align_event;
            figures.SPE.panel_config(p).xdata = {core_data.spikes.(align_event).time_vector};
            figures.SPE.panel_config(p).ydata = {psth_data_spe{i_align}{1}};
            figures.SPE.panel_config(p).c_lims = c_lims;
            p = p+1;
        end
        % Row 2: Heatmaps (cond 2)
        for i_align = 1:numel(align_events)
            align_event = align_events{i_align};
            figures.SPE.panel_config(p).dataType = 'Heatmap';
            figures.SPE.panel_config(p).alignEvent = align_event;
            figures.SPE.panel_config(p).xdata = {core_data.spikes.(align_event).time_vector};
            figures.SPE.panel_config(p).ydata = {psth_data_spe{i_align}{2}};
            figures.SPE.panel_config(p).c_lims = c_lims;
            p = p+1;
        end
        % Row 3: PSTHs
        for i_align = 1:numel(align_events)
            align_event = align_events{i_align};
            figures.SPE.panel_config(p).dataType = 'PSTH';
            figures.SPE.panel_config(p).alignEvent = align_event;
            figures.SPE.panel_config(p).xdata = {core_data.spikes.(align_event).time_vector, core_data.spikes.(align_event).time_vector};
            figures.SPE.panel_config(p).ydata = grand_mean_data_spe{i_align};
            figures.SPE.panel_config(p).colors = cell2mat(color_sets_spe{i_align}');
            figures.SPE.panel_config(p).legendStrings = legend_sets_spe{i_align};
            p = p+1;
        end
        % Row 4: Pupil
        for i_align = 1:numel(align_events)
            align_event = align_events{i_align};
            figures.SPE.panel_config(p).dataType = 'Pupil';
            figures.SPE.panel_config(p).alignEvent = align_event;
            figures.SPE.panel_config(p).xdata = {core_data.pupil.(align_event).time_vector, core_data.pupil.(align_event).time_vector};
            mean_trace1 = mean(core_data.pupil.(align_event).traces(pupil_conditions_spe{i_align}{1}, :), 1, 'omitnan');
            mean_trace2 = mean(core_data.pupil.(align_event).traces(pupil_conditions_spe{i_align}{2}, :), 1, 'omitnan');
            figures.SPE.panel_config(p).ydata = {mean_trace1, mean_trace2};
            figures.SPE.panel_config(p).colors = cell2mat(color_sets_spe{i_align}');
            p = p+1;
        end
    end
end

%% GENERATE PLOTS
giveFeed('Step 7: Generating plots...');

figure_names = fieldnames(figures);
for i_fig = 1:numel(figure_names)
    fig_name = figure_names{i_fig};
    fig_config = figures.(fig_name);

    if isempty(fig_config.panel_config)
        continue;
    end

    fig_h = figure('Position', [100, 100, 1200, 700], ...
        'Color', 'w', 'MenuBar', 'None', 'ToolBar', 'None');
    fig_h.PaperPositionMode = 'auto';

    % --- Layout Calculation ---
    n_cols = 3;
    n_rows = 4; % 2xHeatmap, 1xPSTH, 1xPupil

    % Proportional widths
    time_windows = [diff(core_data.spikes.CUE_ON.window), ...
                    diff(core_data.spikes.outcomeOn.window), ...
                    diff(core_data.spikes.reward.window)];
    proportions = time_windows / sum(time_windows);
    left_margin = 0.07; right_margin = 0.05; h_gap = 0.06;
    plot_area_width = 1 - left_margin - right_margin - (n_cols-1)*h_gap;
    col_widths = plot_area_width * proportions;

    col_lefts = zeros(1, n_cols);
    col_lefts(1) = left_margin;
    for i_col = 2:n_cols
        col_lefts(i_col) = col_lefts(i_col-1) + col_widths(i_col-1) + h_gap;
    end

    % Fixed heights
    row_proportions = [0.15, 0.15, 0.3, 0.3]; % Relative heights
    top_margin = 0.1; bottom_margin = 0.1; v_gap = 0.08;
    plot_area_height = 1 - top_margin - bottom_margin - (n_rows-1)*v_gap;
    row_heights = plot_area_height * row_proportions;

    row_bottoms = zeros(1, n_rows);
    row_bottoms(n_rows) = bottom_margin;
    for i_row = (n_rows-1):-1:1
        row_bottoms(i_row) = row_bottoms(i_row+1) + ...
            row_heights(i_row+1) + v_gap;
    end

    % --- Plotting Loop ---
    for i_panel = 1:numel(fig_config.panel_config)
        panel = fig_config.panel_config(i_panel);

        col = rem(i_panel-1, n_cols) + 1;
        row = floor((i_panel-1) / n_cols) + 1;

        ax_pos = [col_lefts(col), row_bottoms(row), ...
            col_widths(col), row_heights(row)];
        ax = axes('Position', ax_pos);

        ax.Tag = sprintf('%s_%s_Axis', panel.alignEvent, panel.dataType);

        hold on;

        switch panel.dataType
            case 'Heatmap'
                imagesc(ax, panel.xdata{1}, ...
                    1:size(panel.ydata{1}, 1), panel.ydata{1});
                clim(ax, panel.c_lims);
                colormap(ax, flipud(bone(256)));
            case 'PSTH'
                for i_trace = 1:numel(panel.ydata)
                    hb = barStairsFill(panel.xdata{i_trace}, ...
                        zeros(size(panel.ydata{i_trace})), ...
                        panel.ydata{i_trace});
                    set(hb(3), 'Color', panel.colors(i_trace, :));
                    delete(hb(1:2)); % Remove edge and baseline
                end
                if col == n_cols && ~isempty(panel.legendStrings)
                    legend(ax, panel.legendStrings, 'Location', 'best', ...
                        'Interpreter', 'none');
                end
            case 'Pupil'
                for i_trace = 1:numel(panel.ydata)
                    plot(ax, panel.xdata{i_trace}, ...
                        panel.ydata{i_trace}, ...
                        'Color', panel.colors(i_trace, :), ...
                        'LineWidth', 1.5);
                end
        end

        xlim(ax, core_data.spikes.(panel.alignEvent).window);
        grid(ax, 'on');
        xline(ax, 0, 'k--', 'LineWidth', 1);
        set(ax, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 1);
    end

    % --- Post-Plot Formatting ---
    all_axes = findobj(fig_h, 'Type', 'axes');

    % De-clutter tick labels and add labels
    for i_ax = 1:numel(all_axes)
        ax = all_axes(i_ax);
        pos = get(ax, 'Position');

        is_left_col = abs(pos(1) - col_lefts(1)) < 1e-4;

        current_row = 0;
        for i_row = 1:n_rows
            if abs(pos(2) - row_bottoms(i_row)) < 1e-4
                current_row = i_row;
                break;
            end
        end

        if ~is_left_col
            set(ax, 'YTickLabel', []);
        else
            if current_row == 1; ylabel(ax, 'Neurons (Cond 1)'); end
            if current_row == 2; ylabel(ax, 'Neurons (Cond 2)'); end
            if current_row == 3; ylabel(ax, 'Firing Rate (spikes/s)'); end
            if current_row == 4; ylabel(ax, 'Pupil Size (norm.)'); end
        end

        if current_row ~= n_rows
            set(ax, 'XTickLabel', []);
        else
            xlabel(ax, ['Time from ' strrep(ax.Tag, '_', ' ') ' (s)']);
        end
    end

    % Link Y-axes
    psth_axes = findobj(fig_h, 'Type', 'axes', '-regexp', 'Tag', '.*_PSTH_Axis');
    if numel(psth_axes) > 1; linkaxes(psth_axes, 'y'); end

    pupil_axes = findobj(fig_h, 'Type', 'axes', '-regexp', 'Tag', '.*_Pupil_Axis');
    if numel(pupil_axes) > 1; linkaxes(pupil_axes, 'y'); end

    % Add title and save
    sgtitle(fig_h, fig_config.title, 'Interpreter', 'none');
    giveFeed(['Saving ' fig_name ' verification figure...']);
    pdfSave(fig_config.fileName, fig_h.Position(3:4)/72, fig_h);
    giveFeed([fig_name ' figure saved to: ' fig_config.fileName]);
end

%% Finalize and Save Manifest
giveFeed('Step 7: Saving updated manifest...');
writetable(manifest, manifest_path);
giveFeed('Manifest saved. Script finished.');
end

%% Helper function
function [grand_mean_psth, per_neuron_psth] = calculate_mean_psth(rates, ...
    condition_mask)
%calculate_mean_psth Calculates the mean PSTH, excluding silent trials.
%   This function calculates the grand-average PSTH across all neurons and
%   the average PSTH for each individual neuron. It filters out trials
%   where a neuron had a total spike count of zero across the analysis
%   window on a per-neuron basis.
%
%   Args:
%       rates: An [nNeurons x nTrials x nTimeBins] matrix of spike rates.
%       condition_mask: A logical vector indicating which trials to include.
%
%   Returns:
%       grand_mean_psth: A [1 x nTimeBins] vector of the grand-average PSTH.
%       per_neuron_psth: An [nNeurons x nTimeBins] matrix of per-neuron PSTHs.

    n_neurons = size(rates, 1);
    n_time_bins = size(rates, 3);

    per_neuron_psth = zeros(n_neurons, n_time_bins);

    for i_neuron = 1:n_neurons
        % Get the rates for the current neuron and the specified condition
        neuron_rates = squeeze(rates(i_neuron, condition_mask, :));

        % Find trials where the neuron was not silent
        active_trials = sum(neuron_rates, 2, 'omitnan') > 0;

        % Calculate the mean PSTH for the active trials
        if any(active_trials)
            per_neuron_psth(i_neuron, :) = mean(neuron_rates( ...
                active_trials, :), 1, 'omitnan');
        else
            % If the neuron is silent in all trials, its mean PSTH is zero
            per_neuron_psth(i_neuron, :) = zeros(1, n_time_bins);
        end
        per_neuron_psth(i_neuron, :) = mean(neuron_rates, 1, 'omitnan');

    end

    % get rid of rows composed exclusively of 0 or NaN:
    per_neuron_psth(isnan(per_neuron_psth)) = 0;
    per_neuron_psth(all(per_neuron_psth == 0,2), :) = [];

    % Calculate the grand-average PSTH from the per-neuron averages
    grand_mean_psth = mean(per_neuron_psth, 1, 'omitnan');
end
