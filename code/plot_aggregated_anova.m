%% plot_aggregated_anova.m
%
% Generates a summary figure visualizing the proportion of significant
% neurons over time for each ANOVA model term, comparing SC and SNc
% populations.
%
% INPUTS:
%   aggregated_sc_data  - A struct containing aggregated data for SC.
%   aggregated_snc_data - A struct containing aggregated data for SNc.
%
% Author: Jules
% Date: 2025-09-12
%

function plot_aggregated_anova(aggregated_sc_data, aggregated_snc_data)

%% Setup Paths
% Add the 'utils' directory to the path so that helper functions can be
% found.
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Identify ANOVA result fields
% Find all fields in the anova_results struct that start with 'p_'.
% This makes the function robust to different numbers of terms.
all_fields = fieldnames(aggregated_sc_data.anova_results);
anova_fields = all_fields(startsWith(all_fields, 'p_'));
n_terms = numel(anova_fields);

if n_terms == 0
    warning('No ANOVA result fields (starting with ''p_'') found. Cannot plot.');
    return;
end

%% Setup Figure
% Create a multi-panel figure to display the results.
% A 2x3 layout is used for the 6 main model terms.
figure('Position', [100, 100, 1200, 800]);
n_rows = 2;
n_cols = 3;

% Define colors for SC and SNc populations for consistency.
sc_color = [0, 0.4470, 0.7410];  % Blue
snc_color = [0.8500, 0.3250, 0.0980]; % Orange

%% Data Processing and Plotting Loop
h_axes = gobjects(n_terms, 1); % Handles to axes for de-cluttering

for i_term = 1:n_terms
    field_name = anova_fields{i_term};

    % --- Select Subplot ---
    h_axes(i_term) = subplot(n_rows, n_cols, i_term);
    hold on;

    % --- Process and Plot SC Data ---
    p_values_sc = aggregated_sc_data.anova_results.(field_name);
    % Calculate proportion of neurons with a significant effect (p < 0.05)
    % The use of nanmean is critical to handle NaN values from sessions
    % where a given analysis was not applicable (e.g., flicker in main task).
    prop_sig_sc = nanmean(p_values_sc < 0.05, 1);
    n_time_bins = size(prop_sig_sc, 2);
    time_vector = 1:n_time_bins; % Generic time bins

    plot(time_vector, prop_sig_sc, 'Color', sc_color, 'LineWidth', 2);

    % --- Process and Plot SNc Data ---
    p_values_snc = aggregated_snc_data.anova_results.(field_name);
    prop_sig_snc = nanmean(p_values_snc < 0.05, 1);

    plot(time_vector, prop_sig_snc, 'Color', snc_color, 'LineWidth', 2);

    % --- Format Subplot ---
    % Create a clean title from the field name
    clean_title = strrep(field_name, 'p_', '');
    clean_title = strrep(clean_title, '_', ' x ');
    title(clean_title, 'Interpreter', 'none', 'FontWeight', 'normal');

    box off;

    % Set y-axis to be consistent across plots
    ylim([0, 0.5]);
end

%% De-clutter Axes and Add Global Labels
% Adhere to plotting conventions by removing redundant labels.
for i_ax = 1:numel(h_axes)
    % Remove X-tick labels from top row plots
    if i_ax <= n_cols
       set(h_axes(i_ax), 'XTickLabel', []);
    end
    % Remove Y-tick labels from plots not in the first column
    if mod(i_ax, n_cols) ~= 1
        set(h_axes(i_ax), 'YTickLabel', []);
    end
end

% Add global axis labels using a hidden "big" axis
han = axes('Visible','off');
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
ylabel(han, 'Proportion of Neurons (p < 0.05)');
xlabel(han, 'Time Bins');

%% Add Legend
% Create a single, clear legend for the entire figure.
lgd = legend(h_axes(1), {'SC', 'SNc'}, 'Location', 'northeast');
lgd.Box = 'off';

%% Add Overarching Title
% The aggregated data does not specify the alignment event.
% An sgtitle should be added by the calling script to provide this context.
% For example: sgtitle('ANOVA Results Aligned to Outcome');

end
