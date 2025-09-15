%% plot_aggregated_anova.m
%
% Generates a summary figure visualizing the proportion of significant
% neurons over time for each ANOVA model term across different alignment
% events. It compares SC and SNc populations.
%
% This function is designed to work with a nested data structure where
% results are organized first by alignment event, then by p-value type:
% aggregated_data.anova_results.(alignment_event).(p_value_field)
%
% INPUTS:
%   aggregated_data     - A struct containing aggregated data for a brain area.
%   brain_area_name     - A char string (e.g., 'SC' or 'SNc') for labeling.
%
% Author: Jules
% Date: 2025-09-15
%
function plot_aggregated_anova(aggregated_data, brain_area_name)

%% Setup Paths
% Define the project root, ensure the figures directory exists, and add the
% 'utils' directory to the path so that helper functions can be found.
project_root = fullfile(findOneDrive, 'Code', 'tokens-analysis-pipeline');
figures_dir = fullfile(project_root, 'figures');
addpath(fullfile(project_root, 'code', 'utils'));

%% Dynamically Discover Alignment Events and ANOVA Terms
% The function discovers the structure of the analysis from the data itself,
% making it robust to different analysis plans.

% Get alignment events (columns of the plot)
alignment_events = fieldnames(aggregated_data.anova_results);
n_cols = numel(alignment_events);

% Get a unique superset of all p-value fields (rows of the plot)
all_p_fields = {};
for i_event = 1:n_cols
    event_name = alignment_events{i_event};
    field_names = fieldnames(aggregated_data.anova_results.(event_name));
    % Exclude the 'time_vector' field from being treated as a plottable term
    is_p_field = ~strcmp(field_names, 'time_vector');
    all_p_fields = [all_p_fields; field_names(is_p_field)];
end
p_value_fields = unique(all_p_fields);
n_rows = numel(p_value_fields);

if n_rows == 0 || n_cols == 0
    warning('No ANOVA results found in the aggregated data. Cannot plot.');
    return;
end

%% Setup Figure
fig = figure('Position', [100, 100, 1200, 900], 'Color', 'w');
plot_color = [0, 0.4470, 0.7410]; % Blue for SC, adaptable for others
h_axes = gobjects(n_rows, n_cols); % Use a 2D array for axes handles

%% Nested Plotting Loop (Row = ANOVA term, Col = Alignment Event)
for i_row = 1:n_rows
    p_value_name = p_value_fields{i_row};

    for i_col = 1:n_cols
        event_name = alignment_events{i_col};

        % --- Select Subplot ---
        plot_idx = (i_row - 1) * n_cols + i_col;
        h_axes(i_row, i_col) = mySubPlot([n_rows, n_cols, plot_idx]);
        hold on;

        % --- Data Processing and Plotting ---
        % A single time vector is shared for all subplots in a column.
        time_vector = aggregated_data.anova_results.(event_name).time_vector;

        % Check if this specific p-value exists for this event
        if isfield(aggregated_data.anova_results.(event_name), p_value_name)
            % Process and Plot Data
            p_values = aggregated_data.anova_results.(event_name).(p_value_name);
            count_sig = sum(p_values < 0.05, 1, 'omitnan');
            h = barStairsFill(time_vector, zeros(size(count_sig)), count_sig);
            delete(h(2)); % No baseline
            set(h(1), 'FaceColor', plot_color, 'EdgeColor', 'none');
            set(h(3), 'Color', plot_color);
        else
            % If data doesn't exist, display a note on the plot
            text(0.5, 0.5, 'N/A', 'HorizontalAlignment', 'center');
            axis off;
        end
        
        box off;
    end
end

%% De-clutter Axes and Add Labels
for i_row = 1:n_rows
    for i_col = 1:n_cols
        ax = h_axes(i_row, i_col);
        
        % Add Column Titles (Alignment Events) to the top row
        if i_row == 1
            title(ax, strrep(alignment_events{i_col}, '_', ' '), ...
                'FontWeight', 'normal');
        end

        % Add Row Labels (ANOVA Terms) to the first column
        if i_col == 1
            clean_label = strrep(p_value_fields{i_row}, 'p_', '');
            clean_label = strrep(clean_label, '_', ' x ');
            ylabel(ax, clean_label, 'Interpreter', 'none');
        end

        % Remove X-tick labels from all but the bottom row
        if i_row < n_rows
           set(ax, 'XTickLabel', []);
        else
           xlabel_str = sprintf('Time from %s (s)', strrep(event_name, '_', ' '));
           xlabel(ax, xlabel_str);
        end
        
        % Remove Y-tick labels from all but the first column
        if i_col > 1
            set(ax, 'YTickLabel', []);
        end
    end
end

%% Final Figure Touches
% Set common y-limits and styles for all valid axes
all_valid_axes = h_axes(isgraphics(h_axes));
if ~isempty(all_valid_axes)
    [~, yLims] = outerLims(all_valid_axes);
    set(all_valid_axes, 'YLim', yLims, 'TickDir', 'Out', 'Color', 'none', ...
        'XColor', 'k', 'YColor', 'k', 'LineWidth', 1);
end

% Add an overarching title for the entire figure
title_str = sprintf('ANOVA Results for %s: Count of Significant Neurons', ...
    brain_area_name);
sgtitle(title_str, 'FontWeight', 'bold', 'Interpreter', 'none');

% Save figure:
fig_filename = fullfile(figures_dir, ...
    sprintf('aggregated_anova_%s.pdf', brain_area_name));
pdfSave(fig_filename, fig.Position(3:4)/72, fig);

end