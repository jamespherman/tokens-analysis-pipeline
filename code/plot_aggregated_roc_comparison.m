%% plot_aggregated_roc_comparison.m
%
%   Generates a publication-quality summary figure that directly compares
%   SC and SNc population results from the bin-by-bin ROC comparison analysis.
%
% INPUTS:
%   aggregated_sc_data  - A struct containing aggregated data for SC.
%   aggregated_snc_data - A struct containing aggregated data for SNc.
%
% Author: Jules
% Date: 2025-09-12
%

function plot_aggregated_roc_comparison(aggregated_sc_data, aggregated_snc_data)

%% Setup Paths
% Add the 'utils' directory to the path so that helper functions can be
% found.
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Figure and Plotting Setup
fig = figure('Position', [100, 100, 1200, 800], 'Color', 'w');
h_axes = gobjects(2, 3);

% Define the three comparisons to plot
comparisons = struct(...
    'name', {'Dist_at_Cue', 'RPE_at_Outcome', 'RPE_at_Reward'}, ...
    'title', {'Preference: Normal vs. Uniform at Cue Onset', ...
              'Preference: Common vs. Rare High Reward at Outcome', ...
              'Preference: Common vs. Rare High Reward at Reward'}, ...
    'xlabel', {'Time from Cue Onset (s)', 'Time from Outcome Onset (s)', 'Time from Reward (s)'} ...
);

% Define colors for the two populations
sc_color = [0, 0.4470, 0.7410];  % Blue
snc_color = [0.8500, 0.3250, 0.0980]; % Orange
plot_alpha = 0.5; % Transparency for overlapping plots

% we want the right-most axes to be wider than the other two
nCols = [4 4 2];
plotIdx = [1 2 2; 5 6 4];

%% Main Plotting Loop
for i_comp = 1:length(comparisons)
    comp_name = comparisons(i_comp).name;

    % --- Data Validation ---
    if ~isfield(aggregated_sc_data, 'roc_comparison') || ~isfield(aggregated_sc_data.roc_comparison, comp_name)
        warning('plot_aggregated_roc_comparison:no_sc_data', 'No data for %s in SC struct.', comp_name);
        continue;
    end
    if ~isfield(aggregated_snc_data, 'roc_comparison') || ~isfield(aggregated_snc_data.roc_comparison, comp_name)
        warning('plot_aggregated_roc_comparison:no_snc_data', 'No data for %s in SNc struct.', comp_name);
        continue;
    end

    % --- Time Vector Loading ---
    % Load a single session_data.mat file to get the time vectors, which are
    % consistent across all sessions.
    one_drive_path = findOneDrive;
    first_session_id = aggregated_sc_data.session_id{1}; % Use first SC session as template
    session_data_path = fullfile(one_drive_path, 'Neuronal Data Analysis', ...
        first_session_id, [first_session_id '_session_data.mat']);

    if ~exist(session_data_path, 'file')
        error('Could not find session file for %s to get time vectors.', first_session_id);
    end
    temp_data = load(session_data_path, 'session_data');
    time_vector = temp_data.session_data.analysis.roc_comparison.(comp_name).time_vector;

    % --- Data Extraction and Count Calculation ---
    sig_sc = aggregated_sc_data.roc_comparison.(comp_name).sig;
    n_total_sc = size(sig_sc, 1);
    prop_sc_cond2 = sum(sig_sc == 1, 1);
    prop_sc_cond1 = -sum(sig_sc == -1, 1);

    sig_snc = aggregated_snc_data.roc_comparison.(comp_name).sig;
    n_total_snc = size(sig_snc, 1);
    prop_snc_cond2 = sum(sig_snc == 1, 1);
    prop_snc_cond1 = -sum(sig_snc == -1, 1);

    % --- Plotting SC Data (Top Row) ---
    h_axes(1, i_comp) = mySubPlot([2, nCols(i_comp), plotIdx(1,i_comp)]);
    hold on;

    % Plot SC proportions
    h_sc1 = barStairsFill(time_vector, zeros(size(prop_sc_cond2)), prop_sc_cond2);
    delete(h_sc1(2)); % Delete baseline
    set(h_sc1(1), 'FaceColor', sc_color, 'FaceAlpha', plot_alpha, 'EdgeColor', 'none');
    set(h_sc1(3), 'Color', sc_color);

    h_sc2 = barStairsFill(time_vector, zeros(size(prop_sc_cond1)), prop_sc_cond1);
    delete(h_sc2(2)); % Delete baseline
    set(h_sc2(1), 'FaceColor', sc_color, 'FaceAlpha', plot_alpha, 'EdgeColor', 'none');
    set(h_sc2(3), 'Color', sc_color);

    % Formatting for SC plots
    title(h_axes(1, i_comp), comparisons(i_comp).title);
    xlim([time_vector(1), time_vector(end)]);
    line(xlim, [0, 0], 'Color', 'k', 'LineStyle', '--');
    line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--');

    % --- Plotting SNc Data (Bottom Row) ---
    h_axes(2, i_comp) = mySubPlot([2, nCols(i_comp), plotIdx(2,i_comp)]);
    hold on;

    % Plot SNc proportions
    h_snc1 = barStairsFill(time_vector, zeros(size(prop_snc_cond2)), prop_snc_cond2);
    delete(h_snc1(2));
    set(h_snc1(1), 'FaceColor', snc_color, 'FaceAlpha', plot_alpha, 'EdgeColor', 'none');
    set(h_snc1(3), 'Color', snc_color);

    h_snc2 = barStairsFill(time_vector, zeros(size(prop_snc_cond1)), prop_snc_cond1);
    delete(h_snc2(2));
    set(h_snc2(1), 'FaceColor', snc_color, 'FaceAlpha', plot_alpha, 'EdgeColor', 'none');
    set(h_snc2(3), 'Color', snc_color);

    % Formatting for SNc plots
    xlim([time_vector(1), time_vector(end)]);
    line(xlim, [0, 0], 'Color', 'k', 'LineStyle', '--');
    line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--');
end

%% Figure Cleanup and Final Touches
% De-clutter axes per AGENTS.md instructions
set(h_axes(1, :), 'XTickLabel', []); % Remove x-labels from top row
set(h_axes(:, 2:3), 'YTickLabel', []); % Remove y-labels from middle and right columns

% Add axis labels to outer plots
ylabel(h_axes(1, 1), 'Count of Neurons (SC)');
ylabel(h_axes(2, 1), 'Count of Neurons (SNc)');
for i = 1:3
    xlabel(h_axes(2, i), comparisons(i).xlabel);
end

% set axes properties
allAx = findall(fig, 'Type', 'Axes');

% set common y-limits:
[~, yLims] = outerLims(allAx);
set(allAx, 'YLim', yLims, 'TickDir', 'Out');

% Create a single legend for the entire figure
h_sc_patch = patch(NaN, NaN, sc_color, 'FaceAlpha', plot_alpha);
h_snc_patch = patch(NaN, NaN, snc_color, 'FaceAlpha', plot_alpha);
legend([h_sc_patch, h_snc_patch], {'SC', 'SNc'}, ...
    'location', 'best', 'Box', 'off');

sgtitle('Aggregated Population Preference: SC vs. SNc', ...
    'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none');

end
