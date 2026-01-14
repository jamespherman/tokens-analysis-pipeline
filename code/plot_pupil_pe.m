%% plot_pupil_pe.m
%
% Generates a figure visualizing pupil responses to RPE, SPE, and their
% interaction. Layout depends on whether session has AV trials.
%
% AV Sessions (4 columns x 2 rows):
%   Row 1: Mean pupil traces with 95% CI
%     Panel A: RPE effect (4 traces: Rare-Low, Common-Low, Rare-High, Common-High)
%     Panel B: SPE effect (4 traces)
%     Panel C: RPE x SPE for Low rewards (4 traces)
%     Panel D: RPE x SPE for High rewards (4 traces)
%   Row 2: P-value time courses
%     Panel E: 2-way ANOVA p-values for RPE (p_dist, p_mag, p_interaction)
%     Panel F: p_contrast (ranksum: surprising vs certain)
%     Panel G: 3-way ANOVA p-values relevant to Low rewards
%     Panel H: 3-way ANOVA p-values relevant to High rewards
%
% No-AV Sessions (1 column x 2 rows):
%   Panel A: RPE traces (4 conditions)
%   Panel E: RPE p-values (2-way ANOVA)
%
% INPUTS:
%   pupil_pe_results - Struct from analyze_pupil_pe.m (includes is_av_session flag)
%   session_id       - String identifier for the session (for title)
%
% OUTPUT:
%   fig - Handle to the created figure
%
% Author: Claude
% Date: 2026-01-12

function fig = plot_pupil_pe(pupil_pe_results, session_id)

%% Setup Paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);
figures_dir = fullfile(project_root, 'figures');

%% Extract Common Data
time_vec = pupil_pe_results.time_vector;
alpha_thresh = 0.05;
is_av_session = pupil_pe_results.is_av_session;

%% Define Color Palettes

% Panel A: RPE colors (4 conditions: Rare-Low, Common-Low, Rare-High, Common-High)
% Blue for Rare (Normal), Orange for Common (Uniform)
% Light shades for Low, Dark shades for High
rpe_colors = [ ...
    0.40, 0.60, 0.85; ...  % Rare-Low (light blue)
    0.95, 0.60, 0.30; ...  % Common-Low (light orange)
    0.10, 0.25, 0.60; ...  % Rare-High (dark blue)
    0.80, 0.35, 0.10];     % Common-High (dark orange)

% Panel B: SPE colors (categorical)
spe_colors = [ ...
    0.50, 0.50, 0.50; ...  % Gray - No-Flicker (0%)
    0.80, 0.60, 0.20; ...  % Gold - Flicker Omitted (50%)
    0.90, 0.30, 0.10; ...  % Red-orange - Flicker Surprising (50%)
    0.20, 0.70, 0.40];     % Green - Flicker Certain (100%)

% Panel C/D: Interaction colors for Low and High rewards
% Use solid for Expected, dashed for Unexpected
% Blue for Rare, Orange for Common
int_colors_low = [ ...
    0.40, 0.60, 0.85; ...  % Rare-Low, Expected (light blue)
    0.40, 0.60, 0.85; ...  % Rare-Low, Unexpected (light blue, dashed)
    0.95, 0.60, 0.30; ...  % Common-Low, Expected (light orange)
    0.95, 0.60, 0.30];     % Common-Low, Unexpected (light orange, dashed)

int_colors_high = [ ...
    0.10, 0.25, 0.60; ...  % Rare-High, Expected (dark blue)
    0.10, 0.25, 0.60; ...  % Rare-High, Unexpected (dark blue, dashed)
    0.80, 0.35, 0.10; ...  % Common-High, Expected (dark orange)
    0.80, 0.35, 0.10];     % Common-High, Unexpected (dark orange, dashed)

% P-value panel colors
p_color_dist = [0.20, 0.40, 0.70];     % Blue for Distribution
p_color_mag = [0.60, 0.20, 0.60];      % Purple for Magnitude
p_color_int = [0.80, 0.00, 0.60];      % Magenta for Interaction
p_color_spe = [0.20, 0.70, 0.40];      % Green for SPE
p_color_contrast = [0.90, 0.30, 0.10]; % Red-orange for SPE contrast
p_color_dist_spe = [0.00, 0.50, 0.50]; % Teal for Dist x SPE

%% Create Figure with Layout Based on Session Type
if is_av_session
    % AV session: 4 columns x 2 rows
    fig = figure('Position', [50, 50, 1800, 700], 'Color', 'w');
    panel_width = 0.20;
    panel_height = 0.35;
    gap_x = 0.04;
    gap_y = 0.08;
    margin_left = 0.05;
    margin_bottom = 0.10;
    x_positions = margin_left + (0:3) * (panel_width + gap_x);
    y_top = margin_bottom + panel_height + gap_y;
    y_bottom = margin_bottom;
else
    % No-AV session: 1 column x 2 rows (RPE only)
    fig = figure('Position', [50, 50, 500, 700], 'Color', 'w');
    panel_width = 0.75;
    panel_height = 0.35;
    gap_y = 0.08;
    margin_left = 0.15;
    margin_bottom = 0.10;
    x_positions = margin_left;
    y_top = margin_bottom + panel_height + gap_y;
    y_bottom = margin_bottom;
end

%% ========== Panel A (Top): RPE Traces ==========
ax_rpe = axes('Position', [x_positions(1), y_top, panel_width, panel_height]);
hold on;

rpe_means = pupil_pe_results.rpe.means;
rpe_ci_lower = pupil_pe_results.rpe.ci_95(:, :, 1);
rpe_ci_upper = pupil_pe_results.rpe.ci_95(:, :, 2);

% Plot CI shading
for i_cond = 1:4
    x_patch = [time_vec, fliplr(time_vec)];
    y_patch = [rpe_ci_lower(i_cond, :), fliplr(rpe_ci_upper(i_cond, :))];
    valid = ~isnan(y_patch);
    if any(valid)
        fill(x_patch(valid), y_patch(valid), rpe_colors(i_cond, :), ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end
end

% Plot mean lines
h_rpe = gobjects(4, 1);
for i_cond = 1:4
    h_rpe(i_cond) = plot(time_vec, rpe_means(i_cond, :), ...
        'Color', rpe_colors(i_cond, :), 'LineWidth', 2);
end

xline(0, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');

ylabel('Pupil Response (norm.)');
xlim([time_vec(1), time_vec(end)]);
box off;
set(ax_rpe, 'TickDir', 'out', 'XTickLabel', []);

legend(h_rpe, pupil_pe_results.rpe.labels, 'Location', 'northwest', ...
    'FontSize', 7, 'Box', 'off');

tc = pupil_pe_results.rpe.trial_counts;
title_str = sprintf('A. RPE: n=[%d,%d,%d,%d]', tc(1), tc(2), tc(3), tc(4));
title(title_str, 'FontWeight', 'normal');

%% ========== Panel E (Bottom): RPE P-values ==========
ax_rpe_p = axes('Position', [x_positions(1), y_bottom, panel_width, panel_height]);
hold on;

p_dist = pupil_pe_results.rpe.p_dist;
p_mag = pupil_pe_results.rpe.p_mag;
p_dist_mag = pupil_pe_results.rpe.p_dist_mag_interaction;

% Shade significant regions
shade_significant_regions(time_vec, p_dist, alpha_thresh, p_color_dist, 0.15);
shade_significant_regions(time_vec, p_mag, alpha_thresh, p_color_mag, 0.15);
shade_significant_regions(time_vec, p_dist_mag, alpha_thresh, p_color_int, 0.15);

% Plot p-value lines
h_p_rpe(1) = plot(time_vec, p_dist, 'Color', p_color_dist, 'LineWidth', 1.5);
h_p_rpe(2) = plot(time_vec, p_mag, 'Color', p_color_mag, 'LineWidth', 1.5);
h_p_rpe(3) = plot(time_vec, p_dist_mag, 'Color', p_color_int, 'LineWidth', 1.5);

% Reference lines
yline(alpha_thresh, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
xline(0, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');

xlabel('Time from Outcome (s)');
ylabel('p-value');
xlim([time_vec(1), time_vec(end)]);
ylim([0, 1]);
box off;
set(ax_rpe_p, 'TickDir', 'out');

legend(h_p_rpe, {'Dist', 'Mag', 'Dist x Mag'}, 'Location', 'northeast', ...
    'FontSize', 7, 'Box', 'off');

title('E. RPE 2-way ANOVA', 'FontWeight', 'normal');

%% ========== Panels B-H: SPE and Interaction (AV sessions only) ==========
if is_av_session

%% ========== Panel B (Top): SPE Traces ==========
ax_spe = axes('Position', [x_positions(2), y_top, panel_width, panel_height]);
hold on;

spe_means = pupil_pe_results.spe.means;
spe_ci_lower = pupil_pe_results.spe.ci_95(:, :, 1);
spe_ci_upper = pupil_pe_results.spe.ci_95(:, :, 2);
p_contrast = pupil_pe_results.spe.p_contrast;

% Shade significant contrast regions at bottom
shade_significant_regions_below(time_vec, p_contrast, alpha_thresh, ...
    p_color_contrast, 0.3, ax_spe);

% Plot CI shading
for i_cond = 1:4
    x_patch = [time_vec, fliplr(time_vec)];
    y_patch = [spe_ci_lower(i_cond, :), fliplr(spe_ci_upper(i_cond, :))];
    valid = ~isnan(y_patch);
    if any(valid)
        fill(x_patch(valid), y_patch(valid), spe_colors(i_cond, :), ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end
end

% Plot mean lines - thicker for key contrast
h_spe = gobjects(4, 1);
line_widths = [1.5, 1.5, 2.5, 2.5];
for i_cond = 1:4
    h_spe(i_cond) = plot(time_vec, spe_means(i_cond, :), ...
        'Color', spe_colors(i_cond, :), 'LineWidth', line_widths(i_cond));
end

xline(0, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');

xlim([time_vec(1), time_vec(end)]);
box off;
set(ax_spe, 'TickDir', 'out', 'XTickLabel', [], 'YTickLabel', []);

spe_legend = {'0% (No-Flicker)', '50% (Omitted)', ...
    '50% (Surprising)', '100% (Certain)'};
legend(h_spe, spe_legend, 'Location', 'northwest', ...
    'FontSize', 7, 'Box', 'off');

tc = pupil_pe_results.spe.trial_counts;
title_str = sprintf('B. SPE: n=[%d,%d,%d,%d]', tc(1), tc(2), tc(3), tc(4));
title(title_str, 'FontWeight', 'normal');

%% ========== Panel F (Bottom): SPE P-values ==========
ax_spe_p = axes('Position', [x_positions(2), y_bottom, panel_width, panel_height]);
hold on;

% Shade significant regions
shade_significant_regions(time_vec, p_contrast, alpha_thresh, ...
    p_color_contrast, 0.2);

% Plot p-value line
plot(time_vec, p_contrast, 'Color', p_color_contrast, 'LineWidth', 1.5);

% Reference lines
yline(alpha_thresh, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
xline(0, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');

xlabel('Time from Outcome (s)');
xlim([time_vec(1), time_vec(end)]);
ylim([0, 1]);
box off;
set(ax_spe_p, 'TickDir', 'out', 'YTickLabel', []);

title('F. SPE Contrast (ranksum)', 'FontWeight', 'normal');

%% ========== Panel C (Top): RPE x SPE Interaction - Low Rewards ==========
ax_int_low = axes('Position', [x_positions(3), y_top, panel_width, panel_height]);
hold on;

% Extract Low reward conditions (indices 1-4 in interaction results)
% Order: Rare-Low Expected, Rare-Low Unexpected, Common-Low Expected, Common-Low Unexpected
int_means = pupil_pe_results.interaction.means;
int_ci_lower = pupil_pe_results.interaction.ci_95(:, :, 1);
int_ci_upper = pupil_pe_results.interaction.ci_95(:, :, 2);

low_indices = [1, 2, 3, 4];
low_means = int_means(low_indices, :);
low_ci_lower = int_ci_lower(low_indices, :);
low_ci_upper = int_ci_upper(low_indices, :);

% Plot CI shading
for i_cond = 1:4
    x_patch = [time_vec, fliplr(time_vec)];
    y_patch = [low_ci_lower(i_cond, :), fliplr(low_ci_upper(i_cond, :))];
    valid = ~isnan(y_patch);
    if any(valid)
        fill(x_patch(valid), y_patch(valid), int_colors_low(i_cond, :), ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end
end

% Plot mean lines (Solid = Expected, Dashed = Unexpected)
h_int_low = gobjects(4, 1);
line_styles_int = {'-', '--', '-', '--'};
for i_cond = 1:4
    h_int_low(i_cond) = plot(time_vec, low_means(i_cond, :), ...
        'Color', int_colors_low(i_cond, :), 'LineWidth', 2, ...
        'LineStyle', line_styles_int{i_cond});
end

xline(0, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');

xlim([time_vec(1), time_vec(end)]);
box off;
set(ax_int_low, 'TickDir', 'out', 'XTickLabel', [], 'YTickLabel', []);

legend_labels_low = {'Rare-Low Exp', 'Rare-Low Unexp', ...
    'Common-Low Exp', 'Common-Low Unexp'};
legend(h_int_low, legend_labels_low, 'Location', 'northwest', ...
    'FontSize', 6, 'Box', 'off');

tc = pupil_pe_results.interaction.trial_counts;
title_str = sprintf('C. RPE x SPE (Low): n=[%d,%d,%d,%d]', ...
    tc(1), tc(2), tc(3), tc(4));
title(title_str, 'FontWeight', 'normal');

%% ========== Panel G (Bottom): 3-way ANOVA p-values (Low-relevant) ==========
ax_int_low_p = axes('Position', [x_positions(3), y_bottom, panel_width, panel_height]);
hold on;

p_dist_3way = pupil_pe_results.interaction.p_dist;
p_spe_3way = pupil_pe_results.interaction.p_spe;
p_dist_spe_3way = pupil_pe_results.interaction.p_dist_spe;

% Shade significant regions
shade_significant_regions(time_vec, p_dist_3way, alpha_thresh, p_color_dist, 0.15);
shade_significant_regions(time_vec, p_spe_3way, alpha_thresh, p_color_spe, 0.15);
shade_significant_regions(time_vec, p_dist_spe_3way, alpha_thresh, p_color_dist_spe, 0.15);

% Plot p-value lines
h_p_low(1) = plot(time_vec, p_dist_3way, 'Color', p_color_dist, 'LineWidth', 1.5);
h_p_low(2) = plot(time_vec, p_spe_3way, 'Color', p_color_spe, 'LineWidth', 1.5);
h_p_low(3) = plot(time_vec, p_dist_spe_3way, 'Color', p_color_dist_spe, 'LineWidth', 1.5);

% Reference lines
yline(alpha_thresh, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
xline(0, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');

xlabel('Time from Outcome (s)');
xlim([time_vec(1), time_vec(end)]);
ylim([0, 1]);
box off;
set(ax_int_low_p, 'TickDir', 'out', 'YTickLabel', []);

legend(h_p_low, {'Dist', 'SPE', 'Dist x SPE'}, 'Location', 'northeast', ...
    'FontSize', 7, 'Box', 'off');

title('G. 3-way ANOVA (Low)', 'FontWeight', 'normal');

%% ========== Panel D (Top): RPE x SPE Interaction - High Rewards ==========
ax_int_high = axes('Position', [x_positions(4), y_top, panel_width, panel_height]);
hold on;

% Extract High reward conditions (indices 5-8 in interaction results)
% Order: Rare-High Expected, Rare-High Unexpected, Common-High Expected, Common-High Unexpected
high_indices = [5, 6, 7, 8];
high_means = int_means(high_indices, :);
high_ci_lower = int_ci_lower(high_indices, :);
high_ci_upper = int_ci_upper(high_indices, :);

% Plot CI shading
for i_cond = 1:4
    x_patch = [time_vec, fliplr(time_vec)];
    y_patch = [high_ci_lower(i_cond, :), fliplr(high_ci_upper(i_cond, :))];
    valid = ~isnan(y_patch);
    if any(valid)
        fill(x_patch(valid), y_patch(valid), int_colors_high(i_cond, :), ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end
end

% Plot mean lines (Solid = Expected, Dashed = Unexpected)
h_int_high = gobjects(4, 1);
for i_cond = 1:4
    h_int_high(i_cond) = plot(time_vec, high_means(i_cond, :), ...
        'Color', int_colors_high(i_cond, :), 'LineWidth', 2, ...
        'LineStyle', line_styles_int{i_cond});
end

xline(0, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');

xlim([time_vec(1), time_vec(end)]);
box off;
set(ax_int_high, 'TickDir', 'out', 'XTickLabel', [], 'YTickLabel', []);

legend_labels_high = {'Rare-High Exp', 'Rare-High Unexp', ...
    'Common-High Exp', 'Common-High Unexp'};
legend(h_int_high, legend_labels_high, 'Location', 'northwest', ...
    'FontSize', 6, 'Box', 'off');

title_str = sprintf('D. RPE x SPE (High): n=[%d,%d,%d,%d]', ...
    tc(5), tc(6), tc(7), tc(8));
title(title_str, 'FontWeight', 'normal');

%% ========== Panel H (Bottom): 3-way ANOVA p-values (High-relevant) ==========
ax_int_high_p = axes('Position', [x_positions(4), y_bottom, panel_width, panel_height]);
hold on;

% Same p-values from 3-way ANOVA, displayed for symmetry with High panel
% Shade significant regions
shade_significant_regions(time_vec, p_dist_3way, alpha_thresh, p_color_dist, 0.15);
shade_significant_regions(time_vec, p_spe_3way, alpha_thresh, p_color_spe, 0.15);
shade_significant_regions(time_vec, p_dist_spe_3way, alpha_thresh, p_color_dist_spe, 0.15);

% Plot p-value lines
h_p_high(1) = plot(time_vec, p_dist_3way, 'Color', p_color_dist, 'LineWidth', 1.5);
h_p_high(2) = plot(time_vec, p_spe_3way, 'Color', p_color_spe, 'LineWidth', 1.5);
h_p_high(3) = plot(time_vec, p_dist_spe_3way, 'Color', p_color_dist_spe, 'LineWidth', 1.5);

% Reference lines
yline(alpha_thresh, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
xline(0, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');

xlabel('Time from Outcome (s)');
xlim([time_vec(1), time_vec(end)]);
ylim([0, 1]);
box off;
set(ax_int_high_p, 'TickDir', 'out', 'YTickLabel', []);

legend(h_p_high, {'Dist', 'SPE', 'Dist x SPE'}, 'Location', 'northeast', ...
    'FontSize', 7, 'Box', 'off');

title('H. 3-way ANOVA (High)', 'FontWeight', 'normal');

%% Synchronize Y-axes for trace panels (AV sessions)
all_ci_lower = [rpe_ci_lower(:); spe_ci_lower(:); int_ci_lower(:)];
all_ci_upper = [rpe_ci_upper(:); spe_ci_upper(:); int_ci_upper(:)];
y_min = min(all_ci_lower, [], 'omitnan');
y_max = max(all_ci_upper, [], 'omitnan');
y_range = y_max - y_min;
y_lim = [y_min - 0.1 * y_range, y_max + 0.1 * y_range];

set(ax_rpe, 'YLim', y_lim);
set(ax_spe, 'YLim', y_lim);
set(ax_int_low, 'YLim', y_lim);
set(ax_int_high, 'YLim', y_lim);

end  % End of is_av_session block for Panels B-H

%% Set Y-limits for RPE panel (no-AV sessions set here, AV sessions already set above)
if ~is_av_session
    y_min = min(rpe_ci_lower(:), [], 'omitnan');
    y_max = max(rpe_ci_upper(:), [], 'omitnan');
    y_range = y_max - y_min;
    y_lim = [y_min - 0.1 * y_range, y_max + 0.1 * y_range];
    set(ax_rpe, 'YLim', y_lim);
end

%% Add Supertitle
if is_av_session
    title_str = sprintf('Pupil Responses: %s', session_id);
else
    title_str = sprintf('Pupil Responses (RPE only): %s', session_id);
end
sgtitle(title_str, 'Interpreter', 'none', 'FontWeight', 'bold');

%% Save Figure
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
if is_av_session
    fig_filename = fullfile(figures_dir, ...
        sprintf('%s_pupil_pe.pdf', session_id));
else
    fig_filename = fullfile(figures_dir, ...
        sprintf('%s_pupil_rpe.pdf', session_id));
end
pdfSave(fig_filename, fig.Position(3:4)/72, fig);

fprintf('Saved pupil PE figure to: %s\n', fig_filename);

end

%% Helper Function: Shade Significant Regions (for p-value panels)
function shade_significant_regions(time_vec, p_values, alpha, color, ...
    face_alpha)
% SHADE_SIGNIFICANT_REGIONS Adds shaded patches where p < alpha
%
% Draws filled patches from y=0 to y=alpha in p-value panels.

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

        x_patch = [t_start, t_end, t_end, t_start];
        y_patch = [0, 0, alpha, alpha];

        fill(x_patch, y_patch, color, 'FaceAlpha', face_alpha, ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end

%% Helper Function: Shade Significant Regions Below Axes
function shade_significant_regions_below(time_vec, p_values, alpha, ...
    color, face_alpha, ax)
% SHADE_SIGNIFICANT_REGIONS_BELOW Adds shaded bars at bottom of axes
%
% Places significance indicators at the bottom 5% of the y-axis range.

    is_sig = p_values < alpha;
    if ~any(is_sig)
        return;
    end

    yl = get(ax, 'YLim');
    if all(yl == [0, 1])
        y_bar_height = 0.02;
        y_bar_bottom = yl(1);
    else
        y_range = yl(2) - yl(1);
        y_bar_height = 0.03 * y_range;
        y_bar_bottom = yl(1);
    end

    diff_sig = diff([0, is_sig, 0]);
    starts = find(diff_sig == 1);
    ends = find(diff_sig == -1) - 1;

    for i_region = 1:length(starts)
        idx_start = starts(i_region);
        idx_end = ends(i_region);

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

        x_patch = [t_start, t_end, t_end, t_start];
        y_patch = [y_bar_bottom, y_bar_bottom, ...
            y_bar_bottom + y_bar_height, y_bar_bottom + y_bar_height];

        fill(x_patch, y_patch, color, 'FaceAlpha', face_alpha, ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end
