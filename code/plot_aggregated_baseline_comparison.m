%% plot_aggregated_baseline_comparison.m
%
%   Generates a publication-quality summary figure that directly compares
%   SC and SNc population results from the baseline comparison analysis,
%   with each subplot representing a different alignment event.
%
% INPUTS:
%   aggregated_sc_data  - A struct containing aggregated data for SC.
%   aggregated_snc_data - A struct containing aggregated data for SNc.
%
% Author: Jules
% Date: 2025-09-14
%

function plot_aggregated_baseline_comparison(aggregated_sc_data, aggregated_snc_data)

%% Setup Paths
% Define the project root and add the 'utils' directory to the path.
project_root = fullfile(findOneDrive, 'Code', 'tokens-analysis-pipeline');
figures_dir = fullfile(project_root, 'figures');
addpath(fullfile(project_root, 'code', 'utils'));

%% Figure and Plotting Setup
% Get the list of alignment events from the data structure.
if isfield(aggregated_sc_data, 'baseline_comparison')
    event_names = fieldnames(aggregated_sc_data.baseline_comparison);
else
    event_names = {};
end
n_events = length(event_names);

if n_events == 0
    disp('No baseline_comparison data found to plot.');
    return;
end

fig = figure('Position', [100, 100, 400 * n_events, 600], 'Color', 'w');
h_axes = gobjects(1, n_events);

% Define colors for SC and SNc
colors = richColors();
sc_color = colors(1,:);
snc_color = colors(2,:);
pos_alpha = 0.6; % Face alpha for positive modulation
neg_alpha = 0.3; % Face alpha for negative modulation

%% Main Plotting Loop
for i_event = 1:n_events
    event_name = event_names{i_event};

    % --- Data Validation ---
    if ~isfield(aggregated_snc_data.baseline_comparison, event_name)
        warning('plot_aggregated_baseline_comparison:no_snc_data', ...
            'No data for event %s in SNc struct. Skipping.', event_name);
        continue;
    end

    % Data is nested by condition; select one representative condition
    sc_conditions = fieldnames(aggregated_sc_data.baseline_comparison.(event_name));
    snc_conditions = fieldnames(aggregated_snc_data.baseline_comparison.(event_name));

    if isempty(sc_conditions) || isempty(snc_conditions)
        warning('plot_aggregated_baseline_comparison:no_conditions', ...
            'No conditions found for event %s. Skipping.', event_name);
        continue;
    end

    % For simplicity, use the first condition found
    sc_condition_name = sc_conditions{1};
    snc_condition_name = snc_conditions{1};

    sc_event_data = aggregated_sc_data.baseline_comparison.(event_name).(sc_condition_name);
    snc_event_data = aggregated_snc_data.baseline_comparison.(event_name).(snc_condition_name);

    % --- Time Vector and Data Extraction ---
    % The time vector is stored within the condition's data struct
    time_vector = sc_event_data.time_vector;

    sig_sc = sc_event_data.sig;
    n_total_sc = size(sig_sc, 1);
    prop_sc_increase = sum(sig_sc == 1, 1) / n_total_sc;
    prop_sc_decrease = -sum(sig_sc == -1, 1) / n_total_sc; % Negative

    sig_snc = snc_event_data.sig;
    n_total_snc = size(sig_snc, 1);
    prop_snc_increase = sum(sig_snc == 1, 1) / n_total_snc;
    prop_snc_decrease = -sum(sig_snc == -1, 1) / n_total_snc; % Negative

    % --- Plotting (SC and SNc on same axes) ---
    h_axes(1, i_event) = mySubPlot([1, n_events, i_event]);
    hold on;

    % Plot SC proportions
    h_sc_inc = barStairsFill(time_vector, zeros(size(prop_sc_increase)), prop_sc_increase);
    set(h_sc_inc(1), 'FaceColor', sc_color, 'EdgeColor', 'none', 'FaceAlpha', pos_alpha);
    delete(h_sc_inc(2:3));
    h_sc_dec = barStairsFill(time_vector, zeros(size(prop_sc_decrease)), prop_sc_decrease);
    set(h_sc_dec(1), 'FaceColor', sc_color, 'EdgeColor', 'none', 'FaceAlpha', neg_alpha);
    delete(h_sc_dec(2:3));

    % Plot SNc proportions as lines
    plot(time_vector, prop_snc_increase, 'Color', snc_color, 'LineWidth', 2);
    plot(time_vector, prop_snc_decrease, 'Color', snc_color, 'LineWidth', 2);

    % --- Formatting ---
    title_str = sprintf('Alignment: %s', strrep(event_name, '_', ' '));
    title(h_axes(1, i_event), title_str, 'Interpreter', 'none');
    xlim([time_vector(1), time_vector(end)]);
    line(xlim, [0, 0], 'Color', 'k', 'LineStyle', '--');
    line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--');

    xlabel_str = sprintf('Time from %s (s)', strrep(event_name, '_', ' '));
    xlabel(xlabel_str);
end

%% Figure Cleanup and Final Touches
% De-clutter axes
if n_events > 1
    set(h_axes(1, 2:end), 'YTickLabel', []);
end

ylabel(h_axes(1, 1), 'Proportion of Neurons');

% Create a custom legend
h_leg_sc = patch(NaN, NaN, sc_color, 'FaceAlpha', pos_alpha);
h_leg_snc_line = line(NaN, NaN, 'Color', snc_color, 'LineWidth', 2);
legend([h_leg_sc, h_leg_snc_line], {'SC', 'SNc'}, 'Location', 'northeast', 'Box', 'off');

% Set axes properties
all_valid_axes = h_axes(isgraphics(h_axes));
set(all_valid_axes, 'TickDir', 'Out', 'Color', 'none', ...
    'XColor', 'k', 'YColor', 'k', 'LineWidth', 1);

sgtitle('Proportion of Modulated Neurons (vs. Baseline)', 'Interpreter', 'none');

% Save figure:
fig_filename = fullfile(figures_dir, 'aggregated_baseline_comparison.pdf');
pdfSave(fig_filename, fig.Position(3:4)/72, fig);

end
