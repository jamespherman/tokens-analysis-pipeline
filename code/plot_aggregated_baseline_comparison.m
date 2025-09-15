%% plot_aggregated_baseline_comparison.m
%
%   Generates a publication-quality summary figure that directly compares
%   SC and SNc population results from the baseline comparison analysis,
%   with each subplot representing a different alignment event.
%
% INPUTS:
%   aggregated_data  - A struct containing aggregated data for one brain area.
%   brain_area_name  - A char string (e.g., 'SC' or 'SNc') for labeling.
%
% Author: Jules
% Date: 2025-09-15
%

function plot_aggregated_baseline_comparison(aggregated_data, brain_area_name)

%% Setup Paths
% Define the project root and add the 'utils' directory to the path.
project_root = fullfile(findOneDrive, 'Code', 'tokens-analysis-pipeline');
figures_dir = fullfile(project_root, 'figures');
addpath(fullfile(project_root, 'code', 'utils'));

%% Dynamically Discover Events and Conditions
% Discover alignment events (for columns)
if isfield(aggregated_data, 'baseline_comparison')
    alignment_events = fieldnames(aggregated_data.baseline_comparison);
    n_cols = numel(alignment_events);
else
    alignment_events = {};
    n_cols = 0;
end

% Discover baseline conditions (for rows)
if n_cols > 0
    first_event = alignment_events{1};
    baseline_conditions = fieldnames(aggregated_data.baseline_comparison.(first_event));
    n_rows = numel(baseline_conditions);
else
    baseline_conditions = {};
    n_rows = 0;
end

if n_rows == 0 || n_cols == 0
    disp('No baseline comparison data found to plot.');
    return;
end

%% Figure and Plotting Setup
fig = figure('Position', [100, 100, 350 * n_cols, 250 * n_rows], 'Color', 'w');
h_axes = gobjects(n_rows, n_cols);

% Define colors and alpha for plotting
colors = richColors();
plot_color = colors(1,:);
pos_alpha = 0.6; % Face alpha for positive modulation
neg_alpha = 0.3; % Face alpha for negative modulation

%% Main Plotting Loop (Row = Condition, Col = Event)
for i_row = 1:n_rows
    condition_name = baseline_conditions{i_row};

    for i_col = 1:n_cols
        event_name = alignment_events{i_col};

        % --- Select Subplot ---
        plot_idx = (i_row - 1) * n_cols + i_col;
        h_axes(i_row, i_col) = mySubPlot([n_rows, n_cols, plot_idx]);
        hold on;

        % --- Data Validation and Extraction ---
        if ~isfield(aggregated_data.baseline_comparison.(event_name), condition_name)
            text(0.5, 0.5, 'N/A', 'HorizontalAlignment', 'center');
            axis off;
            continue;
        end

        event_data = aggregated_data.baseline_comparison.(event_name).(condition_name);
        time_vector = event_data.time_vector;
        sig_data = event_data.sig;
        n_total = size(sig_data, 1);

        prop_increase = sum(sig_data == 1, 1) / n_total;
        prop_decrease = -sum(sig_data == -1, 1) / n_total; % Negative

        % --- Plotting ---
        h_inc = barStairsFill(time_vector, zeros(size(prop_increase)), prop_increase);
        set(h_inc(1), 'FaceColor', plot_color, 'EdgeColor', 'none', 'FaceAlpha', pos_alpha);
        delete(h_inc(2:3));

        h_dec = barStairsFill(time_vector, zeros(size(prop_decrease)), prop_decrease);
        set(h_dec(1), 'FaceColor', plot_color, 'EdgeColor', 'none', 'FaceAlpha', neg_alpha);
        delete(h_dec(2:3));

        % --- Formatting ---
        xlim([time_vector(1), time_vector(end)]);
        line(xlim, [0, 0], 'Color', 'k', 'LineStyle', '--');
        line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--');
        box off;
    end
end

%% Figure Cleanup and Final Touches
for i_row = 1:n_rows
    for i_col = 1:n_cols
        ax = h_axes(i_row, i_col);

        % Add Column Titles (Alignment Events) to the top row
        if i_row == 1
            title(ax, strrep(alignment_events{i_col}, '_', ' '));
        end

        % Add Row Labels (Baseline Conditions) to the first column
        if i_col == 1
            ylabel(ax, strrep(baseline_conditions{i_row}, '_', ' '));
        end

        % Remove X-tick labels from all but the bottom row
        if i_row < n_rows
           set(ax, 'XTickLabel', []);
        else
           xlabel(ax, sprintf('Time from %s (s)', strrep(alignment_events{i_col}, '_', ' ')));
        end

        % Remove Y-tick labels from all but the first column
        if i_col > 1
            set(ax, 'YTickLabel', []);
        end
    end
end

% Set common y-limits and styles
all_valid_axes = h_axes(isgraphics(h_axes));
if ~isempty(all_valid_axes)
    [~, yLims] = outerLims(all_valid_axes);
    set(all_valid_axes, 'YLim', yLims, 'TickDir', 'Out', 'Color', 'none', ...
        'XColor', 'k', 'YColor', 'k', 'LineWidth', 1);
end

% Add an overarching title
title_str = sprintf('Proportion of Modulated Neurons (vs. Baseline) for %s', brain_area_name);
sgtitle(title_str, 'Interpreter', 'none', 'FontWeight', 'bold');

% Save figure:
fig_filename = fullfile(figures_dir, ...
    sprintf('aggregated_baseline_comparison_%s.pdf', brain_area_name));
pdfSave(fig_filename, fig.Position(3:4)/72, fig);

end
