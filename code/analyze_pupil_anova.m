%% analyze_pupil_anova.m
%
% Performs a time-varying 2x3 ANOVA (Distribution x AV Probability) on
% pupil data aligned to outcome onset. The analysis bins the continuous
% pupil trace into 100ms bins and runs a 2-way ANOVA at each bin to test
% for main effects and their interaction.
%
% Factors:
%   - Distribution (2 levels): Normal vs. Uniform
%   - AV Probability (3 levels): 0%, 50%, 100%
%
% INPUTS:
%   core_data  - Struct containing processed pupil data with field
%                .pupil.outcomeOn containing .traces (n_trials x n_samples)
%                and .time_vector
%   conditions - Struct of logical masks from define_task_conditions.m
%
% OUTPUT:
%   results - Struct containing:
%       .p_dist           - [1 x nBins] p-values for Distribution main
%                           effect
%       .p_av_prob        - [1 x nBins] p-values for AV Probability main
%                           effect
%       .p_interaction    - [1 x nBins] p-values for interaction
%       .time_vector      - [1 x nBins] bin centers in seconds
%       .means            - [6 x nBins] mean pupil per condition
%       .ci_95            - [6 x nBins x 2] bootstrapped 95% CI per
%                           condition, dim 3 is [lower, upper]
%       .trial_counts     - [2 x 3] matrix of trial counts per cell
%       .condition_labels - cell array documenting condition order
%
% Author: Claude
% Date: 2026-01-08

function results = analyze_pupil_anova(core_data, conditions)

%% Setup Paths
% Add the 'utils' directory to the path so that helper functions can be
% found.
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Define Parameters
bin_width_sec = 0.100; % 100 ms bins
n_bootstrap = 1000;    % Number of bootstrap iterations for CI

%% Extract Pupil Data
pupil_traces = core_data.pupil.outcomeOn.traces;
time_vector_raw = core_data.pupil.outcomeOn.time_vector;
[n_trials, n_samples] = size(pupil_traces);

%% Create Factor Vectors
% Distribution factor (2 levels)
is_normal = conditions.is_normal_dist(:);
is_uniform = conditions.is_uniform_dist(:);

% AV Probability factor (3 levels)
% 0%: no flicker, certain no AV
is_av_0 = conditions.is_noflicker_certain(:);
% 50%: flicker surprising (AV occurred when 50% prob) or flicker omitted
% (no AV when 50% prob)
is_av_50 = conditions.is_flicker_surprising(:) | ...
    conditions.is_flicker_omitted(:);
% 100%: flicker certain (AV always occurs)
is_av_100 = conditions.is_flicker_certain(:);

%% Identify Valid Trials
% Only include trials that belong to exactly one level of each factor
valid_dist = (is_normal | is_uniform) & ~(is_normal & is_uniform);
valid_av = (is_av_0 | is_av_50 | is_av_100);

% Ensure no trial belongs to multiple AV levels
av_level_count = double(is_av_0) + double(is_av_50) + double(is_av_100);
valid_av = valid_av & (av_level_count == 1);

% Combined validity mask
valid_trials = valid_dist & valid_av;

%% Build Factor Vectors for Valid Trials
% Distribution factor: 1 = Normal, 2 = Uniform
dist_factor = nan(n_trials, 1);
dist_factor(is_normal) = 1;
dist_factor(is_uniform) = 2;

% AV Probability factor: 1 = 0%, 2 = 50%, 3 = 100%
av_factor = nan(n_trials, 1);
av_factor(is_av_0) = 1;
av_factor(is_av_50) = 2;
av_factor(is_av_100) = 3;

%% Compute Trial Counts Per Cell
% Create trial count matrix [2 x 3] (Dist x AV)
results.trial_counts = zeros(2, 3);
for i_dist = 1:2
    for i_av = 1:3
        results.trial_counts(i_dist, i_av) = sum(valid_trials & ...
            dist_factor == i_dist & av_factor == i_av);
    end
end

%% Bin the Pupil Data
% Compute bin edges based on time vector
time_start = time_vector_raw(1);
time_end = time_vector_raw(end);
bin_edges = time_start:bin_width_sec:time_end;
n_bins = length(bin_edges) - 1;

% Compute bin centers for output
bin_centers = bin_edges(1:end-1) + bin_width_sec / 2;

% Initialize binned data matrix
binned_pupil = nan(n_trials, n_bins);

% Average samples within each bin
for i_bin = 1:n_bins
    bin_start = bin_edges(i_bin);
    bin_end = bin_edges(i_bin + 1);

    % Find samples within this bin
    in_bin = (time_vector_raw >= bin_start) & ...
        (time_vector_raw < bin_end);

    if any(in_bin)
        % Average across samples in bin, handling NaN values
        binned_pupil(:, i_bin) = mean(pupil_traces(:, in_bin), 2, ...
            'omitnan');
    end
end

%% Pre-allocate Results
results.p_dist = nan(1, n_bins);
results.p_av_prob = nan(1, n_bins);
results.p_interaction = nan(1, n_bins);
results.means = nan(6, n_bins);
% ci_95: [6 x nBins x 2], dim 3 is [lower, upper]
results.ci_95 = nan(6, n_bins, 2);
results.time_vector = bin_centers;

% Document the order of conditions in means/ci_95
results.condition_labels = { ...
    'Norm-0%', 'Norm-50%', 'Norm-100%', ...
    'Uni-0%', 'Uni-50%', 'Uni-100%'};

%% Run ANOVA at Each Time Bin
for i_bin = 1:n_bins
    % Extract pupil values for this bin
    pupil_bin = binned_pupil(:, i_bin);

    % Select only valid trials
    valid_idx = valid_trials & ~isnan(pupil_bin);

    % Get data for valid trials
    y = pupil_bin(valid_idx);
    g_dist = dist_factor(valid_idx);
    g_av = av_factor(valid_idx);

    % Check if we have enough data for ANOVA
    % Need at least 2 levels per factor
    if length(unique(g_dist)) < 2 || length(unique(g_av)) < 2
        continue;
    end

    % Check if we have at least some trials in each cell
    if any(results.trial_counts(:) == 0)
        % If any cell is empty, we can still run ANOVA but results may
        % be unreliable. Continue anyway but ANOVA may fail.
    end

    % Run 2-way ANOVA with interaction
    try
        [~, tbl, ~] = anovan(y, {g_dist, g_av}, ...
            'model', 'interaction', ...
            'varnames', {'Distribution', 'AV_Probability'}, ...
            'display', 'off');

        % Extract p-values from the ANOVA table
        % Row 2: Distribution main effect
        % Row 3: AV Probability main effect
        % Row 4: Interaction
        results.p_dist(i_bin) = tbl{2, 7};
        results.p_av_prob(i_bin) = tbl{3, 7};
        results.p_interaction(i_bin) = tbl{4, 7};

    catch ME
        % If ANOVA fails, leave p-values as NaN
        % This can happen with insufficient data
        % Uncomment for debugging:
        % fprintf('ANOVA failed at bin %d: %s\n', i_bin, ME.message);
    end

    % Compute means and bootstrapped 95% CI for each cell
    % Order: Norm-0%, Norm-50%, Norm-100%, Uni-0%, Uni-50%, Uni-100%
    cell_idx = 1;
    for i_dist = 1:2
        for i_av = 1:3
            cell_mask = valid_idx & ...
                dist_factor == i_dist & av_factor == i_av;
            cell_data = pupil_bin(cell_mask);

            % Remove NaN values for bootstrap
            cell_data = cell_data(~isnan(cell_data));
            n_valid = length(cell_data);

            % Compute mean
            if n_valid > 0
                results.means(cell_idx, i_bin) = mean(cell_data);

                % Compute bootstrapped 95% CI
                if n_valid >= 2
                    boot_means = nan(n_bootstrap, 1);
                    for i_boot = 1:n_bootstrap
                        boot_idx = randi(n_valid, n_valid, 1);
                        boot_means(i_boot) = mean(cell_data(boot_idx));
                    end
                    ci_bounds = prctile(boot_means, [2.5, 97.5]);
                    results.ci_95(cell_idx, i_bin, 1) = ci_bounds(1);
                    results.ci_95(cell_idx, i_bin, 2) = ci_bounds(2);
                else
                    % Not enough data for bootstrap, use mean as both
                    % bounds
                    results.ci_95(cell_idx, i_bin, 1) = ...
                        results.means(cell_idx, i_bin);
                    results.ci_95(cell_idx, i_bin, 2) = ...
                        results.means(cell_idx, i_bin);
                end
            end

            cell_idx = cell_idx + 1;
        end
    end
end

%% Print Summary
fprintf('Pupil ANOVA analysis complete.\n');
fprintf('  Total valid trials: %d / %d\n', sum(valid_trials), n_trials);
fprintf('  Trial counts per cell:\n');
fprintf('              AV 0%%    AV 50%%   AV 100%%\n');
fprintf('    Normal:   %4d      %4d      %4d\n', ...
    results.trial_counts(1, 1), results.trial_counts(1, 2), ...
    results.trial_counts(1, 3));
fprintf('    Uniform:  %4d      %4d      %4d\n', ...
    results.trial_counts(2, 1), results.trial_counts(2, 2), ...
    results.trial_counts(2, 3));

end
