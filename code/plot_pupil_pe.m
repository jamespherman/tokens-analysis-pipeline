%% plot_pupil_pe.m
%
% Generates a six-panel figure (2 rows x 3 columns) visualizing pupil
% responses to RPE, SPE, and their interaction.
%
% Top row: Mean pupil traces with 95% CI
%   Panel A: RPE effect (Normal distribution) - 3 traces
%   Panel B: SPE effect - 4 traces
%   Panel C: RPE x AV Interaction - 4 traces
%
% Bottom row: P-value time courses
%   Panel A: p_rpe (1-way ANOVA)
%   Panel B: p_contrast (ranksum: surprising vs certain)
%   Panel C: p_rpe, p_av, p_int (2-way ANOVA)
%
% INPUTS:
%   pupil_pe_results - Struct from analyze_pupil_pe.m
%   session_id       - String identifier for the session (for title)
%
% OUTPUT:
%   fig - Handle to the created figure
%
% Author: Claude
% Date: 2026-01-09

function fig = plot_pupil_pe(pupil_pe_results, session_id)

%% Setup Paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);
figures_dir = fullfile(project_root, 'figures');

%% Extract Common Data
time_vec = pupil_pe_results.time_vector;
alpha_thresh = 0.05;

%% Define Color Palettes

% Panel A: RPE colors (blue gradient: light to dark for low to high)
rpe_colors = [ ...
    0.60, 0.80, 0.95; ...  % Light blue - Rare-Low
    0.30, 0.50, 0.70; ...  % Medium blue - Common
    0.00, 0.20, 0.50];     % Dark blue - Rare-High

% Panel B: SPE colors (categorical)
spe_colors = [ ...
    0.50, 0.50, 0.50; ...  % Gray - No-Flicker (0%)
    0.80, 0.60, 0.20; ...  % Gold - Flicker Omitted (50%)
    0.90, 0.30, 0.10; ...  % Red-orange - Flicker Surprising (50%)
    0.20, 0.70, 0.40];     % Green - Flicker Certain (100%)

% Panel C: Interaction colors (2x2 design)
int_colors = [ ...
    0.40, 0.60, 0.80; ...  % Light blue - Common, No-AV
    0.80, 0.50, 0.40; ...  % Light red - Common, AV
    0.10, 0.30, 0.60; ...  % Dark blue - Rare-High, No-AV
    0.70, 0.20, 0.20];     % Dark red - Rare-High, AV

% P-value panel colors
p_color_rpe = [0.80, 0.40, 0.00];      % Orange-brown for RPE
p_color_av = [0.00, 0.60, 0.50];       % Teal for AV
p_color_int = [0.80, 0.00, 0.60];      % Magenta for Interaction
p_color_contrast = [0.90, 0.30, 0.10]; % Red-orange for SPE contrast

%% Create Figure (2 rows x 3 columns)
fig = figure('Position', [50, 50, 1400, 700], 'Color', 'w');

%% ========== Panel A (Top): RPE Traces ==========
ax_rpe = mySubPlot([2, 3, 1]);
ax_rpe.Position = [0.06, 0.55, 0.28, 0.35];
hold on;

rpe_means = pupil_pe_results.rpe.means;
rpe_ci_lower = pupil_pe_results.rpe.ci_95(:, :, 1);
rpe_ci_upper = pupil_pe_results.rpe.ci_95(:, :, 2);
p_rpe = pupil_pe_results.rpe.p_rpe;

% Shade significant regions at bottom of trace panel
shade_significant_regions_below(time_vec, p_rpe, alpha_thresh, ...
    [0.3, 0.3, 0.3], 0.3, ax_rpe);

% Plot CI shading
for i_cond = 1:3
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
h_rpe = gobjects(3, 1);
for i_cond = 1:3
    h_rpe(i_cond) = plot(time_vec, rpe_means(i_cond, :), ...
        'Color', rpe_colors(i_cond, :), 'LineWidth', 2);
end

xline(0, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');

ylabel('Pupil Response (norm.)');
xlim([time_vec(1), time_vec(end)]);
box off;
set(ax_rpe, 'TickDir', 'out', 'XTickLabel', []);

legend(h_rpe, pupil_pe_results.rpe.labels, 'Location', 'northwest', ...
    'FontSize', 8, 'Box', 'off');

tc = pupil_pe_results.rpe.trial_counts;
title_str = sprintf('A. RPE (Normal Dist): n=[%d, %d, %d]', ...
    tc(1), tc(2), tc(3));
title(title_str, 'FontWeight', 'normal');

%% ========== Panel A (Bottom): RPE P-values ==========
ax_rpe_p = mySubPlot([2, 3, 4]);
ax_rpe_p.Position = [0.06, 0.10, 0.28, 0.35];
hold on;

% Shade significant regions
shade_significant_regions(time_vec, p_rpe, alpha_thresh, ...
    p_color_rpe, 0.2);

% Plot p-value line
plot(time_vec, p_rpe, 'Color', p_color_rpe, 'LineWidth', 1.5);

% Reference lines
yline(alpha_thresh, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
xline(0, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');

xlabel('Time from Outcome (s)');
ylabel('p-value');
xlim([time_vec(1), time_vec(end)]);
ylim([0, 1]);
box off;
set(ax_rpe_p, 'TickDir', 'out');

title('RPE Effect (1-way ANOVA)', 'FontWeight', 'normal');

%% ========== Panel B (Top): SPE Traces ==========
ax_spe = mySubPlot([2, 3, 2]);
ax_spe.Position = [0.39, 0.55, 0.28, 0.35];
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
    'FontSize', 8, 'Box', 'off');

tc = pupil_pe_results.spe.trial_counts;
title_str = sprintf('B. SPE: n=[%d, %d, %d, %d]', ...
    tc(1), tc(2), tc(3), tc(4));
title(title_str, 'FontWeight', 'normal');

%% ========== Panel B (Bottom): SPE P-values ==========
ax_spe_p = mySubPlot([2, 3, 5]);
ax_spe_p.Position = [0.39, 0.10, 0.28, 0.35];
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

title('Surprising vs Certain (ranksum)', 'FontWeight', 'normal');

%% ========== Panel C (Top): Interaction Traces ==========
ax_int = mySubPlot([2, 3, 3]);
ax_int.Position = [0.72, 0.55, 0.26, 0.35];
hold on;

int_means = pupil_pe_results.interaction.means;
int_ci_lower = pupil_pe_results.interaction.ci_95(:, :, 1);
int_ci_upper = pupil_pe_results.interaction.ci_95(:, :, 2);
p_int = pupil_pe_results.interaction.p_int;

% Shade significant interaction regions at bottom
shade_significant_regions_below(time_vec, p_int, alpha_thresh, ...
    p_color_int, 0.3, ax_int);

% Plot CI shading
for i_cond = 1:4
    x_patch = [time_vec, fliplr(time_vec)];
    y_patch = [int_ci_lower(i_cond, :), fliplr(int_ci_upper(i_cond, :))];
    valid = ~isnan(y_patch);
    if any(valid)
        fill(x_patch(valid), y_patch(valid), int_colors(i_cond, :), ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end
end

% Plot mean lines (Solid = No-AV, Dashed = AV)
h_int = gobjects(4, 1);
line_styles = {'-', '--', '-', '--'};
for i_cond = 1:4
    h_int(i_cond) = plot(time_vec, int_means(i_cond, :), ...
        'Color', int_colors(i_cond, :), 'LineWidth', 2, ...
        'LineStyle', line_styles{i_cond});
end

xline(0, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');

xlim([time_vec(1), time_vec(end)]);
box off;
set(ax_int, 'TickDir', 'out', 'XTickLabel', [], 'YTickLabel', []);

legend(h_int, pupil_pe_results.interaction.labels, ...
    'Location', 'northwest', 'FontSize', 8, 'Box', 'off');

tc = pupil_pe_results.interaction.trial_counts;
title_str = sprintf('C. RPE x AV: n=[%d,%d; %d,%d]', ...
    tc(1, 1), tc(1, 2), tc(2, 1), tc(2, 2));
title(title_str, 'FontWeight', 'normal');

%% ========== Panel C (Bottom): Interaction P-values ==========
ax_int_p = mySubPlot([2, 3, 6]);
ax_int_p.Position = [0.72, 0.10, 0.26, 0.35];
hold on;

p_rpe_int = pupil_pe_results.interaction.p_rpe;
p_av_int = pupil_pe_results.interaction.p_av;
p_int_int = pupil_pe_results.interaction.p_int;

% Shade significant regions for each effect
shade_significant_regions(time_vec, p_rpe_int, alpha_thresh, ...
    p_color_rpe, 0.2);
shade_significant_regions(time_vec, p_av_int, alpha_thresh, ...
    p_color_av, 0.2);
shade_significant_regions(time_vec, p_int_int, alpha_thresh, ...
    p_color_int, 0.2);

% Plot p-value lines
h_p(1) = plot(time_vec, p_rpe_int, 'Color', p_color_rpe, ...
    'LineWidth', 1.5);
h_p(2) = plot(time_vec, p_av_int, 'Color', p_color_av, ...
    'LineWidth', 1.5);
h_p(3) = plot(time_vec, p_int_int, 'Color', p_color_int, ...
    'LineWidth', 1.5);

% Reference lines
yline(alpha_thresh, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
xline(0, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');

xlabel('Time from Outcome (s)');
xlim([time_vec(1), time_vec(end)]);
ylim([0, 1]);
box off;
set(ax_int_p, 'TickDir', 'out', 'YTickLabel', []);

legend(h_p, {'RPE', 'AV', 'Interaction'}, 'Location', 'northeast', ...
    'FontSize', 8, 'Box', 'off');

title('2-way ANOVA p-values', 'FontWeight', 'normal');

%% Synchronize Y-axes for trace panels
all_ci_lower = [rpe_ci_lower(:); spe_ci_lower(:); int_ci_lower(:)];
all_ci_upper = [rpe_ci_upper(:); spe_ci_upper(:); int_ci_upper(:)];
y_min = min(all_ci_lower, [], 'omitnan');
y_max = max(all_ci_upper, [], 'omitnan');
y_range = y_max - y_min;
y_lim = [y_min - 0.1 * y_range, y_max + 0.1 * y_range];

set(ax_rpe, 'YLim', y_lim);
set(ax_spe, 'YLim', y_lim);
set(ax_int, 'YLim', y_lim);

%% Add Supertitle
sgtitle(sprintf('Pupil Responses: %s', session_id), ...
    'Interpreter', 'none', 'FontWeight', 'bold');

%% Save Figure
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
fig_filename = fullfile(figures_dir, ...
    sprintf('%s_pupil_pe.pdf', session_id));
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
