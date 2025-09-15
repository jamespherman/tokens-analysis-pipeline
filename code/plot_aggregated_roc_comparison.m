%% plot_aggregated_roc_comparison.m
%
%   Generates a publication-quality summary figure that directly compares
%   SC and SNc population results from the bin-by-bin ROC comparison analysis,
%   organized in a grid by event and comparison type.
%
% INPUTS:
%   aggregated_data  - A struct containing aggregated data for one brain area.
%   brain_area_name  - A char string (e.g., 'SC' or 'SNc') for labeling.
%
% Author: Jules
% Date: 2025-09-15
%

function plot_aggregated_roc_comparison(aggregated_data, brain_area_name)

%% Setup Paths
% Define the project root and add the 'utils' directory to the path.
project_root = fullfile(findOneDrive, 'Code', 'tokens-analysis-pipeline');
figures_dir = fullfile(project_root, 'figures');
addpath(fullfile(project_root, 'code', 'utils'));

%% Grid and Figure Setup
% Get event names for columns from the data structure.
event_names = fieldnames(aggregated_data.roc_comparison);
n_cols = length(event_names);

% Get all unique comparison names for rows.
all_comp_names = {};
for i_event = 1:n_cols
    event_name = event_names{i_event};
    comps_for_event = fieldnames(aggregated_data.roc_comparison.(event_name));
    all_comp_names = [all_comp_names; comps_for_event];
end
comp_names = unique(all_comp_names, 'stable');
n_rows = length(comp_names);

% Create figure and axes handles for a grid.
fig = figure('Position', [100, 100, 350 * n_cols, 300 * n_rows], 'Color', 'w');
h_axes = gobjects(n_rows, n_cols);

% Define colors for the two preferences
colors = richColors();
posColor = colors(7,:);
negColor = colors(10,:);

%% Main Plotting Loop (Events as columns, Comparisons as rows)
for i_row = 1:n_rows
    comp_name = comp_names{i_row};

    for i_col = 1:n_cols
        event_name = event_names{i_col};

        % --- Select Subplot ---
        plot_idx = (i_row - 1) * n_cols + i_col;
        h_axes(i_row, i_col) = mySubPlot([n_rows, n_cols, plot_idx]);
        hold on;

        % --- Data Validation and Extraction ---
        if ~isfield(aggregated_data.roc_comparison.(event_name), comp_name)
            text(0.5, 0.5, 'N/A', 'HorizontalAlignment', 'center');
            axis off;
            continue;
        end

        comp_data = aggregated_data.roc_comparison.(event_name).(comp_name);

        if ~isfield(comp_data, 'sig')
             warning('plot_aggregated_roc_comparison:missing_data', ...
                'Missing sig data for %s/%s.', event_name, comp_name);
            continue;
        end

        time_vector = comp_data.time_vector;
        sig_data = comp_data.sig;
        count_cond2 = sum(sig_data == 1, 1);
        count_cond1 = -sum(sig_data == -1, 1);

        % --- Plotting Data ---
        h1 = barStairsFill(time_vector, zeros(size(count_cond2)), count_cond2);
        delete(h1(2)); set(h1(1), 'FaceColor', posColor, 'EdgeColor', 'none'); set(h1(3), 'Color', posColor);
        h2 = barStairsFill(time_vector, zeros(size(count_cond1)), count_cond1);
        delete(h2(2)); set(h2(1), 'FaceColor', negColor, 'EdgeColor', 'none'); set(h2(3), 'Color', negColor);
        xlim([time_vector(1), time_vector(end)]);
        line(xlim, [0, 0], 'Color', 'k', 'LineStyle', '--');
        line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--');
        box off;
    end
end

%% Figure Cleanup and Final Touches
% Add titles to columns and labels to rows/columns
for i_col = 1:n_cols
    % Add column titles (event names) to the top-most plot of each column
    if isgraphics(h_axes(1, i_col))
        title(h_axes(1, i_col), strrep(event_names{i_col}, '_', ' '));
    end
    % Add x-axis labels to the bottom-most plot of each column
    if isgraphics(h_axes(n_rows, i_col))
        xlabel(h_axes(n_rows, i_col), sprintf('Time from %s (s)', strrep(event_names{i_col}, '_', ' ')));
    end
end

for i_row = 1:n_rows
    comp_title = get_comp_title(comp_names{i_row});
    % Add y-axis labels to the left-most plots of each row
    if isgraphics(h_axes(i_row, 1))
        ylabel(h_axes(i_row, 1), comp_title);
    end
end

% De-clutter axes per AGENTS.md instructions
set(h_axes(1:n_rows-1, :), 'XTickLabel', []);
if n_cols > 1
    set(h_axes(:, 2:end), 'YTickLabel', []);
end

% Set common y-limits and styles
all_valid_axes = h_axes(isgraphics(h_axes));
if ~isempty(all_valid_axes)
    [~, yLims] = outerLims(all_valid_axes);
    set(all_valid_axes, 'YLim', yLims, 'TickDir', 'Out', 'Color', 'none', ...
        'XColor', 'k', 'YColor', 'k', 'LineWidth', 1);
end

sgtitle(sprintf('Aggregated ROC Comparison for %s', brain_area_name), 'Interpreter', 'none');

% Save figure
fig_filename = fullfile(figures_dir, sprintf('aggregated_roc_comparison_%s.pdf', brain_area_name));
pdfSave(fig_filename, fig.Position(3:4)/72, fig);

end

%% Local Helper Function
function title_str = get_comp_title(comp_name)
    % Generates a descriptive title for a given comparison name.
    switch comp_name
        case 'Dist'
            title_str = 'Norm vs. Unif';
        case 'RPE'
            title_str = 'Common vs. Rare';
        case 'SPE'
            title_str = 'Surprise vs. No Surprise';
        otherwise
            title_str = strrep(comp_name, '_', ' ');
    end
end
