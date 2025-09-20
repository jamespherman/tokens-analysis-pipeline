%% plot_aggregated_svm.m
%
%   Generates a publication-quality summary figure showing time-resolved
%   SVM classification accuracy. It plots individual session traces and the
%   mean trace for a 3x3 grid of comparisons.
%
% INPUTS:
%   aggregated_data  - A struct containing aggregated data for one brain area.
%   brain_area_name  - A char string (e.g., 'SC' or 'SNc') for labeling.
%
% Author: Jules
% Date: 2025-09-15
%

function plot_aggregated_svm(aggregated_data, brain_area_name)

%% Setup Paths
project_root = fullfile(findOneDrive, 'Code', 'tokens-analysis-pipeline');
figures_dir = fullfile(project_root, 'figures');
addpath(fullfile(project_root, 'code', 'utils'));

%% Grid and Figure Setup
event_names = fieldnames(aggregated_data.svm_results);
n_cols = length(event_names);
comp_names = {'Dist'; 'RPE'; 'SPE'};
n_rows = length(comp_names);

fig = figure('Position', [100, 100, 350 * n_cols, 300 * n_rows], 'Color', 'w');
h_axes = gobjects(n_rows, n_cols);

light_gray = [0.8, 0.8, 0.8];
dark_gray = [0.3, 0.3, 0.3];

%% Main Plotting Loop
for i_row = 1:n_rows
    row_type = comp_names{i_row};
    for i_col = 1:n_cols
        event_name = event_names{i_col};
        plot_idx = (i_row - 1) * n_cols + i_col;
        h_axes(i_row, i_col) = mySubPlot([n_rows, n_cols, plot_idx]);
        hold on;

        available_comps = fieldnames(aggregated_data.svm_results.(event_name));
        actual_comp_name_cell = available_comps(startsWith(available_comps, row_type));

        if isempty(actual_comp_name_cell)
            text(0.5, 0.5, 'N/A', 'HorizontalAlignment', 'center');
            axis off;
            continue;
        end

        actual_comp_name = actual_comp_name_cell{1};
        comp_data = aggregated_data.svm_results.(event_name).(actual_comp_name);

        if ~isfield(comp_data, 'accuracy')
            warning('plot_aggregated_svm:missing_data', ...
                'Missing accuracy data for %s/%s.', event_name, actual_comp_name);
            continue;
        end

        time_vector = comp_data.time_vector;
        accuracy_matrix = comp_data.accuracy; % [n_sessions x n_time_bins]

        % Plot individual session traces
        for i_session = 1:size(accuracy_matrix, 1)
            h_indiv = barStairsFill(time_vector, ...
                zeros(size(time_vector)), accuracy_matrix(i_session, :));
            delete(h_indiv(1:2)); % Delete bar and baseline
            set(h_indiv(3), 'Color', light_gray, 'LineWidth', 0.5);
        end

        % Plot mean trace
        mean_accuracy = mean(accuracy_matrix, 1, 'omitnan');
        h_mean = barStairsFill(time_vector, ...
            zeros(size(time_vector)), mean_accuracy);
        delete(h_mean(1:2)); % Delete bar and baseline
        set(h_mean(3), 'Color', dark_gray, 'LineWidth', 2);

        xlim([time_vector(1), time_vector(end)]);
        ylim([0.3 1.0]); % Set y-axis limits for accuracy
        line(xlim, [0.5, 0.5], 'Color', 'k', 'LineStyle', '--');
        line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--');
        box off;
    end
end

%% Figure Cleanup and Final Touches
for i_col = 1:n_cols
    if isgraphics(h_axes(1, i_col))
        title(h_axes(1, i_col), strrep(event_names{i_col}, '_', ' '));
    end
    if isgraphics(h_axes(n_rows, i_col))
        xlabel(h_axes(n_rows, i_col), sprintf('Time from %s (s)', strrep(event_names{i_col}, '_', ' ')));
    end
end

for i_row = 1:n_rows
    if isgraphics(h_axes(i_row, 1))
        ylabel(h_axes(i_row, 1), 'Classification Accuracy');
    end
end

set(h_axes(1:n_rows-1, :), 'XTickLabel', []);
if n_cols > 1
    set(h_axes(:, 2:end), 'YTickLabel', []);
end

all_valid_axes = h_axes(isgraphics(h_axes));
if ~isempty(all_valid_axes)
    set(all_valid_axes, 'TickDir', 'Out', 'Color', 'none', ...
        'XColor', 'k', 'YColor', 'k', 'LineWidth', 1);
end

sgtitle(sprintf('Aggregated SVM Classification for %s', brain_area_name), 'Interpreter', 'none');

fig_filename = fullfile(figures_dir, sprintf('aggregated_svm_classification_%s.pdf', brain_area_name));
pdfSave(fig_filename, fig.Position(3:4)/72, fig);

end
