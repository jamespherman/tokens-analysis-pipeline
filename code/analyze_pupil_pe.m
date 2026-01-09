%% analyze_pupil_pe.m
%
% Performs time-varying analyses of pupil responses to Reward Prediction
% Error (RPE), Sensory Prediction Error (SPE), and their interaction.
%
% Panel A - RPE (Normal distribution only):
%   3 conditions: rare-low, common, rare-high
%   1-way ANOVA at each time bin
%
% Panel B - SPE:
%   4 conditions: noflicker_certain, flicker_omitted, flicker_surprising,
%                 flicker_certain
%   Pairwise contrast: flicker_surprising vs flicker_certain (ranksum)
%
% Panel C - RPE x AV Interaction (Normal distribution only):
%   4 conditions from 2x2 crossing of RPE (common vs rare-high) with
%   AV presence (No-AV vs AV)
%   2-way ANOVA at each time bin
%
% INPUTS:
%   core_data  - Struct containing processed pupil data with field
%                .pupil.outcomeOn containing .traces and .time_vector
%   conditions - Struct of logical masks from define_task_conditions.m
%
% OUTPUT:
%   results - Struct containing:
%       Panel A (RPE):
%           .rpe.p_rpe         - [1 x nBins] p-values for RPE effect
%           .rpe.means         - [3 x nBins] means per condition
%           .rpe.ci_95         - [3 x nBins x 2] bootstrapped 95% CI
%           .rpe.trial_counts  - [1 x 3] trial counts
%           .rpe.labels        - condition labels
%       Panel B (SPE):
%           .spe.p_contrast    - [1 x nBins] p-values for surprising vs
%                                certain (ranksum)
%           .spe.means         - [4 x nBins] means per condition
%           .spe.ci_95         - [4 x nBins x 2] bootstrapped 95% CI
%           .spe.trial_counts  - [1 x 4] trial counts
%           .spe.labels        - condition labels
%       Panel C (Interaction):
%           .interaction.p_rpe - [1 x nBins] p-values for RPE main effect
%           .interaction.p_av  - [1 x nBins] p-values for AV main effect
%           .interaction.p_int - [1 x nBins] p-values for interaction
%           .interaction.means - [4 x nBins] means per condition
%           .interaction.ci_95 - [4 x nBins x 2] bootstrapped 95% CI
%           .interaction.trial_counts - [2 x 2] trial counts
%           .interaction.labels - condition labels
%       .time_vector - [1 x nBins] bin centers in seconds
%
% Author: Claude
% Date: 2026-01-09

function results = analyze_pupil_pe(core_data, conditions)

%% Setup Paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Define Parameters
bin_width_sec = 0.100;
n_bootstrap = 1000;

%% Extract Pupil Data
pupil_traces = core_data.pupil.outcomeOn.traces;
time_vector_raw = core_data.pupil.outcomeOn.time_vector;
[n_trials, ~] = size(pupil_traces);

%% Bin the Pupil Data
time_start = time_vector_raw(1);
time_end = time_vector_raw(end);
bin_edges = time_start:bin_width_sec:time_end;
n_bins = length(bin_edges) - 1;
bin_centers = bin_edges(1:end-1) + bin_width_sec / 2;

binned_pupil = nan(n_trials, n_bins);
for i_bin = 1:n_bins
    bin_start = bin_edges(i_bin);
    bin_end = bin_edges(i_bin + 1);
    in_bin = (time_vector_raw >= bin_start) & (time_vector_raw < bin_end);
    if any(in_bin)
        binned_pupil(:, i_bin) = mean(pupil_traces(:, in_bin), 2, ...
            'omitnan');
    end
end

results.time_vector = bin_centers;

%% ========== Panel A: RPE Analysis (Normal Distribution Only) ==========

% Define RPE conditions (3 levels)
rpe_masks = { ...
    conditions.is_norm_rare_low(:), ...
    conditions.is_norm_common(:), ...
    conditions.is_norm_rare_high(:)};
rpe_labels = {'Rare-Low', 'Common', 'Rare-High'};

% Compute trial counts
results.rpe.trial_counts = cellfun(@sum, rpe_masks);
results.rpe.labels = rpe_labels;

% Pre-allocate
results.rpe.p_rpe = nan(1, n_bins);
results.rpe.means = nan(3, n_bins);
results.rpe.ci_95 = nan(3, n_bins, 2);

% Build RPE factor for ANOVA
rpe_factor = nan(n_trials, 1);
for i_cond = 1:3
    rpe_factor(rpe_masks{i_cond}) = i_cond;
end
valid_rpe = ~isnan(rpe_factor);

% Run 1-way ANOVA at each bin
for i_bin = 1:n_bins
    pupil_bin = binned_pupil(:, i_bin);
    valid_idx = valid_rpe & ~isnan(pupil_bin);

    y = pupil_bin(valid_idx);
    g = rpe_factor(valid_idx);

    if length(unique(g)) >= 2
        try
            [~, tbl, ~] = anovan(y, {g}, 'display', 'off');
            results.rpe.p_rpe(i_bin) = tbl{2, 7};
        catch
            % Leave as NaN
        end
    end

    % Compute means and CI for each condition
    for i_cond = 1:3
        cell_mask = rpe_masks{i_cond} & ~isnan(pupil_bin);
        cell_data = pupil_bin(cell_mask);
        cell_data = cell_data(~isnan(cell_data));
        n_valid = length(cell_data);

        if n_valid > 0
            results.rpe.means(i_cond, i_bin) = mean(cell_data);
            if n_valid >= 2
                boot_means = nan(n_bootstrap, 1);
                for i_boot = 1:n_bootstrap
                    boot_idx = randi(n_valid, n_valid, 1);
                    boot_means(i_boot) = mean(cell_data(boot_idx));
                end
                ci_bounds = prctile(boot_means, [2.5, 97.5]);
                results.rpe.ci_95(i_cond, i_bin, 1) = ci_bounds(1);
                results.rpe.ci_95(i_cond, i_bin, 2) = ci_bounds(2);
            else
                results.rpe.ci_95(i_cond, i_bin, :) = ...
                    results.rpe.means(i_cond, i_bin);
            end
        end
    end
end

%% ========== Panel B: SPE Analysis ==========

% Define SPE conditions (4 levels)
spe_masks = { ...
    conditions.is_noflicker_certain(:), ...
    conditions.is_flicker_omitted(:), ...
    conditions.is_flicker_surprising(:), ...
    conditions.is_flicker_certain(:)};
spe_labels = {'No-Flicker (0%)', 'Flicker Omitted (50%)', ...
    'Flicker Surprising (50%)', 'Flicker Certain (100%)'};

% Compute trial counts
results.spe.trial_counts = cellfun(@sum, spe_masks);
results.spe.labels = spe_labels;

% Pre-allocate
results.spe.p_contrast = nan(1, n_bins);
results.spe.means = nan(4, n_bins);
results.spe.ci_95 = nan(4, n_bins, 2);

% Key contrast: flicker_surprising vs flicker_certain
mask_surprising = spe_masks{3};
mask_certain = spe_masks{4};

% Run ranksum at each bin for the key contrast
for i_bin = 1:n_bins
    pupil_bin = binned_pupil(:, i_bin);

    % Key contrast (ranksum test)
    data_surp = pupil_bin(mask_surprising & ~isnan(pupil_bin));
    data_cert = pupil_bin(mask_certain & ~isnan(pupil_bin));

    if length(data_surp) >= 2 && length(data_cert) >= 2
        try
            results.spe.p_contrast(i_bin) = ranksum(data_surp, data_cert);
        catch
            % Leave as NaN
        end
    end

    % Compute means and CI for each condition
    for i_cond = 1:4
        cell_mask = spe_masks{i_cond} & ~isnan(pupil_bin);
        cell_data = pupil_bin(cell_mask);
        cell_data = cell_data(~isnan(cell_data));
        n_valid = length(cell_data);

        if n_valid > 0
            results.spe.means(i_cond, i_bin) = mean(cell_data);
            if n_valid >= 2
                boot_means = nan(n_bootstrap, 1);
                for i_boot = 1:n_bootstrap
                    boot_idx = randi(n_valid, n_valid, 1);
                    boot_means(i_boot) = mean(cell_data(boot_idx));
                end
                ci_bounds = prctile(boot_means, [2.5, 97.5]);
                results.spe.ci_95(i_cond, i_bin, 1) = ci_bounds(1);
                results.spe.ci_95(i_cond, i_bin, 2) = ci_bounds(2);
            else
                results.spe.ci_95(i_cond, i_bin, :) = ...
                    results.spe.means(i_cond, i_bin);
            end
        end
    end
end

%% ========== Panel C: RPE x AV Interaction (Normal Dist Only) ==========

% AV presence factor:
%   No-AV: noflicker_certain OR flicker_omitted
%   AV: flicker_certain OR flicker_surprising
is_no_av = conditions.is_noflicker_certain(:) | ...
    conditions.is_flicker_omitted(:);
is_av = conditions.is_flicker_certain(:) | ...
    conditions.is_flicker_surprising(:);

% RPE factor for interaction (2 levels: common vs rare-high)
is_common = conditions.is_norm_common(:);
is_rare_high = conditions.is_norm_rare_high(:);

% 4 conditions for 2x2 design:
% 1: Common + No-AV
% 2: Common + AV
% 3: Rare-High + No-AV
% 4: Rare-High + AV
int_masks = { ...
    is_common & is_no_av, ...
    is_common & is_av, ...
    is_rare_high & is_no_av, ...
    is_rare_high & is_av};
int_labels = {'Common, No-AV', 'Common, AV', ...
    'Rare-High, No-AV', 'Rare-High, AV'};

% Compute trial counts (2x2 matrix)
results.interaction.trial_counts = [ ...
    sum(int_masks{1}), sum(int_masks{2}); ...
    sum(int_masks{3}), sum(int_masks{4})];
results.interaction.labels = int_labels;

% Pre-allocate
results.interaction.p_rpe = nan(1, n_bins);
results.interaction.p_av = nan(1, n_bins);
results.interaction.p_int = nan(1, n_bins);
results.interaction.means = nan(4, n_bins);
results.interaction.ci_95 = nan(4, n_bins, 2);

% Build factors for 2-way ANOVA
% RPE factor: 1=common, 2=rare-high
% AV factor: 1=no-AV, 2=AV
rpe_int_factor = nan(n_trials, 1);
rpe_int_factor(is_common) = 1;
rpe_int_factor(is_rare_high) = 2;

av_int_factor = nan(n_trials, 1);
av_int_factor(is_no_av) = 1;
av_int_factor(is_av) = 2;

valid_int = ~isnan(rpe_int_factor) & ~isnan(av_int_factor);

% Run 2-way ANOVA at each bin
for i_bin = 1:n_bins
    pupil_bin = binned_pupil(:, i_bin);
    valid_idx = valid_int & ~isnan(pupil_bin);

    y = pupil_bin(valid_idx);
    g_rpe = rpe_int_factor(valid_idx);
    g_av = av_int_factor(valid_idx);

    if length(unique(g_rpe)) >= 2 && length(unique(g_av)) >= 2
        try
            [~, tbl, ~] = anovan(y, {g_rpe, g_av}, ...
                'model', 'interaction', ...
                'varnames', {'RPE', 'AV'}, ...
                'display', 'off');
            results.interaction.p_rpe(i_bin) = tbl{2, 7};
            results.interaction.p_av(i_bin) = tbl{3, 7};
            results.interaction.p_int(i_bin) = tbl{4, 7};
        catch
            % Leave as NaN
        end
    end

    % Compute means and CI for each condition
    for i_cond = 1:4
        cell_mask = int_masks{i_cond} & ~isnan(pupil_bin);
        cell_data = pupil_bin(cell_mask);
        cell_data = cell_data(~isnan(cell_data));
        n_valid = length(cell_data);

        if n_valid > 0
            results.interaction.means(i_cond, i_bin) = mean(cell_data);
            if n_valid >= 2
                boot_means = nan(n_bootstrap, 1);
                for i_boot = 1:n_bootstrap
                    boot_idx = randi(n_valid, n_valid, 1);
                    boot_means(i_boot) = mean(cell_data(boot_idx));
                end
                ci_bounds = prctile(boot_means, [2.5, 97.5]);
                results.interaction.ci_95(i_cond, i_bin, 1) = ci_bounds(1);
                results.interaction.ci_95(i_cond, i_bin, 2) = ci_bounds(2);
            else
                results.interaction.ci_95(i_cond, i_bin, :) = ...
                    results.interaction.means(i_cond, i_bin);
            end
        end
    end
end

%% Print Summary
fprintf('Pupil PE analysis complete.\n');
fprintf('\n--- Panel A: RPE (Normal Distribution) ---\n');
fprintf('  Trial counts: Rare-Low=%d, Common=%d, Rare-High=%d\n', ...
    results.rpe.trial_counts(1), results.rpe.trial_counts(2), ...
    results.rpe.trial_counts(3));

fprintf('\n--- Panel B: SPE ---\n');
fprintf('  Trial counts: No-Flicker=%d, Omitted=%d, Surprising=%d, ', ...
    results.spe.trial_counts(1), results.spe.trial_counts(2), ...
    results.spe.trial_counts(3));
fprintf('Certain=%d\n', results.spe.trial_counts(4));

fprintf('\n--- Panel C: RPE x AV Interaction ---\n');
fprintf('  Trial counts:\n');
fprintf('                No-AV    AV\n');
fprintf('    Common:     %4d   %4d\n', ...
    results.interaction.trial_counts(1, 1), ...
    results.interaction.trial_counts(1, 2));
fprintf('    Rare-High:  %4d   %4d\n', ...
    results.interaction.trial_counts(2, 1), ...
    results.interaction.trial_counts(2, 2));

end
