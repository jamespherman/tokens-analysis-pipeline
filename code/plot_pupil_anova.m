%% plot_pupil_anova.m
%
% Generates a two-panel figure visualizing the results of the pupil ANOVA
% analysis (Distribution x AV Probability).
%
% Panel A (top): Mean pupil traces for all 6 conditions with shaded 95% CI.
%   - Colors distinguish the 3 AV probability levels (0%, 50%, 100%)
%   - Line styles distinguish distributions (solid=Normal, dashed=Uniform)
%
% Panel B (bottom): P-value time course for Distribution, AV Probability,
%   and their interaction. Regions where p < 0.05 are shaded.
%
% INPUTS:
%   pupil_anova_results - Struct from analyze_pupil_anova.m containing:
%                         .means, .ci_95, .time_vector, .p_dist,
%                         .p_av_prob, .p_interaction, .condition_labels
%   session_id          - String identifier for the session (for title)
%
% OUTPUT:
%   fig - Handle to the created figure
%
% Author: Claude
% Date: 2026-01-08

function fig = plot_pupil_anova(pupil_anova_results, session_id)

%% Setup Paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);
figures_dir = fullfile(project_root, 'figures');

%% Extract Data
time_vec = pupil_anova_results.time_vector;
means = pupil_anova_results.means;
ci_lower = pupil_anova_results.ci_95(:, :, 1);
ci_upper = pupil_anova_results.ci_95(:, :, 2);
p_dist = pupil_anova_results.p_dist;
p_av_prob = pupil_anova_results.p_av_prob;
p_interaction = pupil_anova_results.p_interaction;

%% Define Colors and Styles
% Colors for AV probability levels (0%, 50%, 100%)
% Using a colorblind-friendly palette
av_colors = [ ...
    0.00, 0.45, 0.70; ...  % Blue for 0%
    0.90, 0.60, 0.00; ...  % Orange for 50%
    0.00, 0.62, 0.45];     % Teal for 100%

% Line styles for distribution
line_style_normal = '-';
line_style_uniform = '--';

% P-value colors
p_colors = [ ...
    0.80, 0.40, 0.00; ...  % Orange-brown for Distribution
    0.00, 0.60, 0.50; ...  % Teal for AV Probability
    0.80, 0.00, 0.60];     % Magenta for Interaction

%% Create Figure
fig = figure('Position', [100, 100, 900, 700], 'Color', 'w');

%% Panel A: Mean Pupil Traces with 95% CI
ax_traces = mySubPlot([2, 1, 1]);
% Adjust position to add vertical whitespace between panels
% Format: [left, bottom, width, height]
ax_traces.Position = [0.10, 0.55, 0.85, 0.35];
hold on;

% Condition order: Norm-0%, Norm-50%, Norm-100%, Uni-0%, Uni-50%, Uni-100%
% Plot in order: first shaded CI regions (back), then lines (front)

% Plot shaded CI regions for all conditions
for i_cond = 1:6
    % Determine which AV level (1, 2, or 3)
    av_level = mod(i_cond - 1, 3) + 1;
    color = av_colors(av_level, :);

    % Create patch for CI
    x_patch = [time_vec, fliplr(time_vec)];
    y_patch = [ci_lower(i_cond, :), fliplr(ci_upper(i_cond, :))];

    % Remove NaN segments for clean patches
    valid = ~isnan(y_patch);
    if any(valid)
        fill(x_patch(valid), y_patch(valid), color, ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end
end

% Plot mean lines
h_lines = gobjects(6, 1);
for i_cond = 1:6
    % Determine distribution (1-3 = Normal, 4-6 = Uniform)
    is_uniform = i_cond > 3;
    if is_uniform
        line_style = line_style_uniform;
    else
        line_style = line_style_normal;
    end

    % Determine which AV level (1, 2, or 3)
    av_level = mod(i_cond - 1, 3) + 1;
    color = av_colors(av_level, :);

    h_lines(i_cond) = plot(time_vec, means(i_cond, :), ...
        'Color', color, 'LineStyle', line_style, 'LineWidth', 1.5);
end

% Add vertical line at t=0
xline(0, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');

% Configure axes
ylabel('Pupil Response (norm.)');
xlim([time_vec(1), time_vec(end)]);
box off;
set(ax_traces, 'TickDir', 'out');

% Create legend with custom entries
% Group by AV probability with line style indicator
legend_entries = { ...
    'AV 0% (Normal)', 'AV 50% (Normal)', 'AV 100% (Normal)', ...
    'AV 0% (Uniform)', 'AV 50% (Uniform)', 'AV 100% (Uniform)'};
legend(h_lines, legend_entries, 'Location', 'northeast', ...
    'FontSize', 8, 'Box', 'off');

title('A. Pupil Response by Condition', 'FontWeight', 'normal');

%% Panel B: P-value Time Course
ax_pvals = mySubPlot([2, 1, 2]);
% Adjust position to add vertical whitespace between panels
ax_pvals.Position = [0.10, 0.10, 0.85, 0.35];
hold on;

% Shade significant regions (p < 0.05)
alpha_thresh = 0.05;

% Shade for Distribution
shade_significant_regions(time_vec, p_dist, alpha_thresh, ...
    p_colors(1, :), 0.2);

% Shade for AV Probability
shade_significant_regions(time_vec, p_av_prob, alpha_thresh, ...
    p_colors(2, :), 0.2);

% Shade for Interaction
shade_significant_regions(time_vec, p_interaction, alpha_thresh, ...
    p_colors(3, :), 0.2);

% Plot p-value lines
h_p(1) = plot(time_vec, p_dist, 'Color', p_colors(1, :), ...
    'LineWidth', 1.5);
h_p(2) = plot(time_vec, p_av_prob, 'Color', p_colors(2, :), ...
    'LineWidth', 1.5);
h_p(3) = plot(time_vec, p_interaction, 'Color', p_colors(3, :), ...
    'LineWidth', 1.5);

% Add horizontal line at p = 0.05
yline(alpha_thresh, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

% Add vertical line at t=0
xline(0, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');

% Configure axes
xlabel('Time from Outcome Onset (s)');
ylabel('p-value');
xlim([time_vec(1), time_vec(end)]);
ylim([0, 1]);
box off;
set(ax_pvals, 'TickDir', 'out');

% Legend
legend(h_p, {'Distribution', 'AV Probability', 'Interaction'}, ...
    'Location', 'northeast', 'FontSize', 8, 'Box', 'off');

title('B. ANOVA p-values Over Time', 'FontWeight', 'normal');

%% Add Supertitle
% Include trial counts in subtitle
trial_counts = pupil_anova_results.trial_counts;
subtitle_str = sprintf(['Trials: Norm(0%%=%d, 50%%=%d, 100%%=%d), ' ...
    'Uni(0%%=%d, 50%%=%d, 100%%=%d)'], ...
    trial_counts(1, 1), trial_counts(1, 2), trial_counts(1, 3), ...
    trial_counts(2, 1), trial_counts(2, 2), trial_counts(2, 3));

sgtitle({sprintf('Pupil ANOVA: %s', session_id), subtitle_str}, ...
    'Interpreter', 'none', 'FontWeight', 'bold');

%% Save Figure
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
fig_filename = fullfile(figures_dir, ...
    sprintf('%s_pupil_anova.pdf', session_id));
pdfSave(fig_filename, fig.Position(3:4)/72, fig);

fprintf('Saved pupil ANOVA figure to: %s\n', fig_filename);

end

%% Helper Function: Shade Significant Regions
function shade_significant_regions(time_vec, p_values, alpha, color, ...
    face_alpha)
% SHADE_SIGNIFICANT_REGIONS Adds shaded patches where p < alpha
%
% Finds contiguous regions where p < alpha and draws filled patches
% from y=0 to y=alpha.

    % Find where p < alpha
    is_sig = p_values < alpha;

    if ~any(is_sig)
        return;
    end

    % Find contiguous regions
    diff_sig = diff([0, is_sig, 0]);
    starts = find(diff_sig == 1);
    ends = find(diff_sig == -1) - 1;

    for i_region = 1:length(starts)
        idx_start = starts(i_region);
        idx_end = ends(i_region);

        % Get time bounds (extend slightly for visual clarity)
        if idx_start > 1
            t_start = (time_vec(idx_start) + time_vec(idx_start - 1)) / 2;
        else
            t_start = time_vec(idx_start);
        end

        if idx_end < length(time_vec)
            t_end = (time_vec(idx_end) + time_vec(idx_end + 1)) / 2;
        else
            t_end = time_vec(idx_end);
        end

        % Draw patch from y=0 to y=alpha
        x_patch = [t_start, t_end, t_end, t_start];
        y_patch = [0, 0, alpha, alpha];

        fill(x_patch, y_patch, color, 'FaceAlpha', face_alpha, ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end
