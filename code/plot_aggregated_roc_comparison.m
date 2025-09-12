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
figure('Position', [100, 100, 1200, 400]);

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

%% Main Plotting Loop
for i_comp = 1:length(comparisons)
    comp_name = comparisons(i_comp).name;

    % --- Subplot Setup ---
    mySubPlot(1, 3, i_comp);
    hold on;

    % --- Data Extraction ---
    % Note: This function assumes the aggregation script has been fixed to
    % include 'time_vector' in the aggregated data structures.
    % As of 2025-09-12, it is a known issue that it does not.
    if ~isfield(aggregated_sc_data, 'roc_comparison') || ~isfield(aggregated_sc_data.roc_comparison, comp_name)
        warning('plot_aggregated_roc_comparison:no_sc_data', 'No data for %s in SC struct.', comp_name);
        continue;
    end
     if ~isfield(aggregated_snc_data, 'roc_comparison') || ~isfield(aggregated_snc_data.roc_comparison, comp_name)
        warning('plot_aggregated_roc_comparison:no_snc_data', 'No data for %s in SNc struct.', comp_name);
        continue;
    end

    sig_sc = aggregated_sc_data.roc_comparison.(comp_name).sig;
    sig_snc = aggregated_snc_data.roc_comparison.(comp_name).sig;

    % HACK: The time vector is not currently aggregated. We load a data file
    % from the first session in the aggregation list to get it. This should be
    % fixed in `aggregate_analysis_results.m`.
    if isfield(aggregated_sc_data, 'session_id') && ~isempty(aggregated_sc_data.session_id)
        first_session_id = aggregated_sc_data.session_id{1};
        hack_path = fullfile('data', 'processed', first_session_id, 'analysis_results.mat');
        if exist(hack_path, 'file')
            temp_data = load(hack_path, 'analysis_results');
            time_vector = temp_data.analysis_results.roc_comparison.(comp_name).time_vector;
        else
            error('Could not find analysis file for session %s to get time_vector.', first_session_id);
        end
    else
        error('Cannot determine time_vector because aggregated_sc_data is empty or missing session_id.');
    end


    % --- Proportion Calculation ---
    n_total_sc = size(sig_sc, 1);
    prop_sc_cond2 = sum(sig_sc == 1, 1) / n_total_sc;
    prop_sc_cond1 = -sum(sig_sc == -1, 1) / n_total_sc;

    n_total_snc = size(sig_snc, 1);
    prop_snc_cond2 = sum(sig_snc == 1, 1) / n_total_snc;
    prop_snc_cond1 = -sum(sig_snc == -1, 1) / n_total_snc;

    % --- Plotting ---
    % Plot SC data
    barStairsFill(time_vector, prop_sc_cond2, 'FaceColor', sc_color, 'EdgeColor', 'none', 'FaceAlpha', plot_alpha);
    barStairsFill(time_vector, prop_sc_cond1, 'FaceColor', sc_color, 'EdgeColor', 'none', 'FaceAlpha', plot_alpha);

    % Plot SNc data
    barStairsFill(time_vector, prop_snc_cond2, 'FaceColor', snc_color, 'EdgeColor', 'none', 'FaceAlpha', plot_alpha);
    barStairsFill(time_vector, prop_snc_cond1, 'FaceColor', snc_color, 'EdgeColor', 'none', 'FaceAlpha', plot_alpha);

    % --- Formatting ---
    title(comparisons(i_comp).title);
    xlabel(comparisons(i_comp).xlabel);
    xlim([time_vector(1), time_vector(end)]);
    line(xlim, [0, 0], 'Color', 'k', 'LineStyle', '--');
    line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--');

    % De-clutter axes per AGENTS.md instructions
    if i_comp == 1
        ylabel('Proportion of Neurons');
    else
        set(gca, 'YTickLabel', []);
    end
end

%% Legend and Final Touches
% Create a single legend for the entire figure
h_sc = patch(NaN, NaN, sc_color, 'FaceAlpha', plot_alpha);
h_snc = patch(NaN, NaN, snc_color, 'FaceAlpha', plot_alpha);
legend([h_sc, h_snc], {'SC', 'SNc'}, 'Location', 'northwest', 'Box', 'off');

sgtitle('Aggregated Population Preference: SC vs. SNc', 'FontSize', 16, 'FontWeight', 'bold');

end
