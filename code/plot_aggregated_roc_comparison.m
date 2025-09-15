%% plot_aggregated_roc_comparison.m
%
%   Generates a publication-quality summary figure that directly compares
%   SC and SNc population results from the bin-by-bin ROC comparison analysis,
%   organized in a grid by event and comparison type.
%
% INPUTS:
%   aggregated_sc_data  - A struct containing aggregated data for SC.
%   aggregated_snc_data - A struct containing aggregated data for SNc.
%
% Author: Jules
% Date: 2025-09-14
%

function plot_aggregated_roc_comparison(aggregated_sc_data, aggregated_snc_data)

%% Setup Paths
% Define the project root and add the 'utils' directory to the path.
project_root = fullfile(findOneDrive, 'Code', 'tokens-analysis-pipeline');
figures_dir = fullfile(project_root, 'figures');
addpath(fullfile(project_root, 'code', 'utils'));

%% Grid and Figure Setup
% Get event names for columns from the SC data structure.
event_names = fieldnames(aggregated_sc_data.roc_comparison);
n_events = length(event_names);

% Get all unique comparison names for rows.
all_comp_names = {};
for i_event = 1:n_events
    event_name = event_names{i_event};
    comps_for_event = fieldnames(aggregated_sc_data.roc_comparison.(event_name));
    all_comp_names = [all_comp_names; comps_for_event];
end
comp_names = unique(all_comp_names, 'stable');
n_comps = length(comp_names);

% Create figure and axes handles for a grid.
% Each comparison gets two rows (one for SC, one for SNc).
fig = figure('Position', [100, 100, 350 * n_events, 300 * n_comps], 'Color', 'w');
h_axes = gobjects(2 * n_comps, n_events);

% Define colors for the two preferences
colors = richColors;
posColor = colors(7,:);
negColor = colors(10,:);

%% Main Plotting Loop (Events as columns, Comparisons as rows)
for i_event = 1:n_events
    event_name = event_names{i_event};

    for i_comp = 1:n_comps
        comp_name = comp_names{i_comp};

        % --- Data Validation and Extraction ---
        if ~isfield(aggregated_sc_data.roc_comparison, event_name) || ...
           ~isfield(aggregated_sc_data.roc_comparison.(event_name), comp_name) || ...
           ~isfield(aggregated_snc_data.roc_comparison, event_name) || ...
           ~isfield(aggregated_snc_data.roc_comparison.(event_name), comp_name)
            continue; % Skip if data doesn't exist for this combination
        end

        sc_comp_data = aggregated_sc_data.roc_comparison.(event_name).(comp_name);
        snc_comp_data = aggregated_snc_data.roc_comparison.(event_name).(comp_name);

        if ~isfield(sc_comp_data, 'sig') || ~isfield(snc_comp_data, 'sig')
            warning('plot_aggregated_roc_comparison:missing_data', ...
                'Missing sig data for %s/%s.', event_name, comp_name);
            continue;
        end

        % Correctly access the time vector from within the roc_comparison struct
        time_vector = sc_comp_data.time_vectors.sig;
        sig_sc = sc_comp_data.sig;
        count_sc_cond2 = sum(sig_sc == 1, 1);
        count_sc_cond1 = -sum(sig_sc == -1, 1);

        sig_snc = snc_comp_data.sig;
        count_snc_cond2 = sum(sig_snc == 1, 1);
        count_snc_cond1 = -sum(sig_snc == -1, 1);

        % --- Subplot Index Calculation ---
        sc_row_idx = (i_comp - 1) * 2 + 1;
        snc_row_idx = sc_row_idx + 1;
        sc_subplot_idx = (sc_row_idx - 1) * n_events + i_event;
        snc_subplot_idx = (snc_row_idx - 1) * n_events + i_event;

        % --- Plotting SC Data ---
        h_axes(sc_row_idx, i_event) = mySubPlot([2 * n_comps, n_events, sc_subplot_idx]);
        hold on;
        h_sc1 = barStairsFill(time_vector, zeros(size(count_sc_cond2)), count_sc_cond2);
        delete(h_sc1(2)); set(h_sc1(1), 'FaceColor', posColor, 'EdgeColor', 'none'); set(h_sc1(3), 'Color', posColor);
        h_sc2 = barStairsFill(time_vector, zeros(size(count_sc_cond1)), count_sc_cond1);
        delete(h_sc2(2)); set(h_sc2(1), 'FaceColor', negColor, 'EdgeColor', 'none'); set(h_sc2(3), 'Color', negColor);
        xlim([time_vector(1), time_vector(end)]);
        line(xlim, [0, 0], 'Color', 'k', 'LineStyle', '--');
        line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--');

        % --- Plotting SNc Data ---
        h_axes(snc_row_idx, i_event) = mySubPlot([2 * n_comps, n_events, snc_subplot_idx]);
        hold on;
        h_snc1 = barStairsFill(time_vector, zeros(size(count_snc_cond2)), count_snc_cond2);
        delete(h_snc1(2)); set(h_snc1(1), 'FaceColor', posColor, 'EdgeColor', 'none'); set(h_snc1(3), 'Color', posColor);
        h_snc2 = barStairsFill(time_vector, zeros(size(count_snc_cond1)), count_snc_cond1);
        delete(h_snc2(2)); set(h_snc2(1), 'FaceColor', negColor, 'EdgeColor', 'none'); set(h_snc2(3), 'Color', negColor);
        xlim([time_vector(1), time_vector(end)]);
        line(xlim, [0, 0], 'Color', 'k', 'LineStyle', '--');
        line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--');
    end
end

%% Figure Cleanup and Final Touches
% Add titles to columns and labels to rows/columns
for i_event = 1:n_events
    % Add column titles (event names) to the top-most plot of each column
    title(h_axes(1, i_event), strrep(event_names{i_event}, '_', ' '));
    % Add x-axis labels to the bottom-most plot of each column
    xlabel(h_axes(end, i_event), sprintf('Time from %s (s)', strrep(event_names{i_event}, '_', ' ')));
end

for i_comp = 1:n_comps
    sc_row_idx = (i_comp - 1) * 2 + 1;
    snc_row_idx = sc_row_idx + 1;
    comp_title = get_comp_title(comp_names{i_comp});
    % Add y-axis labels to the left-most plots of each row
    ylabel(h_axes(sc_row_idx, 1), sprintf('SC: %s', comp_title));
    ylabel(h_axes(snc_row_idx, 1), sprintf('SNc: %s', comp_title));
end

% De-clutter axes per AGENTS.md instructions
set(h_axes(1:end-1, :), 'XTickLabel', []);
if n_events > 1
    set(h_axes(:, 2:end), 'YTickLabel', []);
end

% Set common y-limits and styles
all_valid_axes = h_axes(isgraphics(h_axes));
if ~isempty(all_valid_axes)
    [~, yLims] = outerLims(all_valid_axes);
    set(all_valid_axes, 'YLim', yLims, 'TickDir', 'Out', 'Color', 'none', ...
        'XColor', 'k', 'YColor', 'k', 'LineWidth', 1);
end

sgtitle('Aggregated ROC Comparison: SC vs. SNc Population Preference', 'Interpreter', 'none');

% Save figure
fig_filename = fullfile(figures_dir, 'aggregated_roc_comparison.pdf');
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
