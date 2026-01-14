%% analyze_pupil_pe.m
%
% Performs time-varying analyses of pupil responses to Reward Prediction
% Error (RPE), Sensory Prediction Error (SPE), and their interaction.
%
% Panel A - RPE (Magnitude-Matched Cross-Distribution):
%   4 conditions: Rare-Low, Common-Low, Rare-High, Common-High
%   2-way ANOVA at each time bin: Distribution (Rare vs Common) x Magnitude (Low vs High)
%
% Panel B - SPE:
%   4 conditions: noflicker_certain, flicker_omitted, flicker_surprising,
%                 flicker_certain
%   Pairwise contrast: flicker_surprising vs flicker_certain (ranksum)
%
% Panels C/D - RPE x SPE Interaction:
%   8 conditions crossing RPE (Rare-Low, Common-Low, Rare-High, Common-High)
%   with SPE (Unexpected AV = is_flicker_surprising, Expected AV = is_flicker_certain)
%   3-way ANOVA: Distribution x SPE x Magnitude
%
% INPUTS:
%   core_data     - Struct containing processed pupil data with field
%                   .pupil.outcomeOn containing .traces and .time_vector
%   conditions    - Struct of logical masks from define_task_conditions.m
%   is_av_session - Boolean flag indicating whether session has AV trials.
%                   If true, SPE and interaction analyses are computed.
%                   If false, only RPE analysis is computed.
%
% OUTPUT:
%   results - Struct containing:
%       Panel A (RPE):
%           .rpe.p_dist           - [1 x nBins] p-values for Distribution effect
%           .rpe.p_mag            - [1 x nBins] p-values for Magnitude effect
%           .rpe.p_dist_mag_interaction - [1 x nBins] p-values for interaction
%           .rpe.means            - [4 x nBins] means per condition
%           .rpe.ci_95            - [4 x nBins x 2] bootstrapped 95% CI
%           .rpe.trial_counts     - [1 x 4] trial counts
%           .rpe.labels           - condition labels
%       Panel B (SPE):
%           .spe.p_contrast       - [1 x nBins] p-values for surprising vs
%                                   certain (ranksum)
%           .spe.means            - [4 x nBins] means per condition
%           .spe.ci_95            - [4 x nBins x 2] bootstrapped 95% CI
%           .spe.trial_counts     - [1 x 4] trial counts
%           .spe.labels           - condition labels
%       Panels C/D (Interaction):
%           .interaction.p_dist       - [1 x nBins] Distribution main effect
%           .interaction.p_spe        - [1 x nBins] SPE main effect
%           .interaction.p_mag        - [1 x nBins] Magnitude main effect
%           .interaction.p_dist_spe   - [1 x nBins] Distribution x SPE
%           .interaction.p_dist_mag   - [1 x nBins] Distribution x Magnitude
%           .interaction.p_spe_mag    - [1 x nBins] SPE x Magnitude
%           .interaction.p_three_way  - [1 x nBins] 3-way interaction
%           .interaction.means        - [8 x nBins] means per condition
%           .interaction.ci_95        - [8 x nBins x 2] bootstrapped 95% CI
%           .interaction.trial_counts - [1 x 8] trial counts
%           .interaction.labels       - condition labels
%       .time_vector - [1 x nBins] bin centers in seconds
%
% Author: Claude
% Date: 2026-01-12

function results = analyze_pupil_pe(core_data, conditions, is_av_session)

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
results.is_av_session = is_av_session;

%% ========== Panel A: RPE Analysis (Magnitude-Matched) ==========

% Define RPE conditions (4 levels: 2 Distributions x 2 Magnitudes)
% Order: Rare-Low, Common-Low, Rare-High, Common-High
rpe_masks = { ...
    conditions.is_rare_low(:), ...
    conditions.is_common_low(:), ...
    conditions.is_rare_high(:), ...
    conditions.is_common_high(:)};
rpe_labels = {'Rare-Low', 'Common-Low', 'Rare-High', 'Common-High'};

% Compute trial counts
results.rpe.trial_counts = cellfun(@sum, rpe_masks);
results.rpe.labels = rpe_labels;

% Pre-allocate
results.rpe.p_dist = nan(1, n_bins);
results.rpe.p_mag = nan(1, n_bins);
results.rpe.p_dist_mag_interaction = nan(1, n_bins);
results.rpe.means = nan(4, n_bins);
results.rpe.ci_95 = nan(4, n_bins, 2);

% Build factors for 2-way ANOVA
% Distribution: 1=Rare (Normal), 2=Common (Uniform)
% Magnitude: 1=Low, 2=High
dist_factor = nan(n_trials, 1);
mag_factor = nan(n_trials, 1);

% Rare-Low: Dist=1, Mag=1
dist_factor(rpe_masks{1}) = 1;
mag_factor(rpe_masks{1}) = 1;

% Common-Low: Dist=2, Mag=1
dist_factor(rpe_masks{2}) = 2;
mag_factor(rpe_masks{2}) = 1;

% Rare-High: Dist=1, Mag=2
dist_factor(rpe_masks{3}) = 1;
mag_factor(rpe_masks{3}) = 2;

% Common-High: Dist=2, Mag=2
dist_factor(rpe_masks{4}) = 2;
mag_factor(rpe_masks{4}) = 2;

valid_rpe = ~isnan(dist_factor) & ~isnan(mag_factor);

% Run 2-way ANOVA at each bin
for i_bin = 1:n_bins
    pupil_bin = binned_pupil(:, i_bin);
    valid_idx = valid_rpe & ~isnan(pupil_bin);

    y = pupil_bin(valid_idx);
    g_dist = dist_factor(valid_idx);
    g_mag = mag_factor(valid_idx);

    if length(unique(g_dist)) >= 2 && length(unique(g_mag)) >= 2
        try
            [~, tbl, ~] = anovan(y, {g_dist, g_mag}, ...
                'model', 'interaction', ...
                'varnames', {'Dist', 'Mag'}, ...
                'display', 'off');
            results.rpe.p_dist(i_bin) = tbl{2, 7};
            results.rpe.p_mag(i_bin) = tbl{3, 7};
            results.rpe.p_dist_mag_interaction(i_bin) = tbl{4, 7};
        catch
            % Leave as NaN
        end
    end

    % Compute means and CI for each condition
    for i_cond = 1:4
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

if is_av_session
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
else
    % No-AV session: set SPE results to empty/NaN
    results.spe.trial_counts = nan(1, 4);
    results.spe.labels = {};
    results.spe.p_contrast = nan(1, n_bins);
    results.spe.means = nan(4, n_bins);
    results.spe.ci_95 = nan(4, n_bins, 2);
end

%% ========== Panels C/D: RPE x SPE Interaction (3-way ANOVA) ==========

if is_av_session
    % SPE factor:
    %   Unexpected AV: flicker_surprising
    %   Expected AV: flicker_certain
    is_unexpected_spe = conditions.is_flicker_surprising(:);
    is_expected_spe = conditions.is_flicker_certain(:);

    % 8 conditions crossing RPE (4 levels) with SPE (2 levels):
    % Order matches: RPE condition x SPE condition
    % 1: Rare-Low x Expected
    % 2: Rare-Low x Unexpected
    % 3: Common-Low x Expected
    % 4: Common-Low x Unexpected
    % 5: Rare-High x Expected
    % 6: Rare-High x Unexpected
    % 7: Common-High x Expected
    % 8: Common-High x Unexpected
    int_masks = { ...
        conditions.is_rare_low(:) & is_expected_spe, ...
        conditions.is_rare_low(:) & is_unexpected_spe, ...
        conditions.is_common_low(:) & is_expected_spe, ...
        conditions.is_common_low(:) & is_unexpected_spe, ...
        conditions.is_rare_high(:) & is_expected_spe, ...
        conditions.is_rare_high(:) & is_unexpected_spe, ...
        conditions.is_common_high(:) & is_expected_spe, ...
        conditions.is_common_high(:) & is_unexpected_spe};
    int_labels = { ...
        'Rare-Low, Expected', 'Rare-Low, Unexpected', ...
        'Common-Low, Expected', 'Common-Low, Unexpected', ...
        'Rare-High, Expected', 'Rare-High, Unexpected', ...
        'Common-High, Expected', 'Common-High, Unexpected'};

    % Compute trial counts
    results.interaction.trial_counts = cellfun(@sum, int_masks);
    results.interaction.labels = int_labels;

    % Pre-allocate p-values for 3-way ANOVA
    results.interaction.p_dist = nan(1, n_bins);
    results.interaction.p_spe = nan(1, n_bins);
    results.interaction.p_mag = nan(1, n_bins);
    results.interaction.p_dist_spe = nan(1, n_bins);
    results.interaction.p_dist_mag = nan(1, n_bins);
    results.interaction.p_spe_mag = nan(1, n_bins);
    results.interaction.p_three_way = nan(1, n_bins);
    results.interaction.means = nan(8, n_bins);
    results.interaction.ci_95 = nan(8, n_bins, 2);

    % Build factors for 3-way ANOVA
    % Distribution: 1=Rare (Normal), 2=Common (Uniform)
    % SPE: 1=Expected, 2=Unexpected
    % Magnitude: 1=Low, 2=High
    dist_int_factor = nan(n_trials, 1);
    spe_int_factor = nan(n_trials, 1);
    mag_int_factor = nan(n_trials, 1);

    for i_cond = 1:8
        mask = int_masks{i_cond};
        % Determine Distribution factor (Rare=1, Common=2)
        if i_cond <= 2  % Rare-Low
            dist_int_factor(mask) = 1;
            mag_int_factor(mask) = 1;
        elseif i_cond <= 4  % Common-Low
            dist_int_factor(mask) = 2;
            mag_int_factor(mask) = 1;
        elseif i_cond <= 6  % Rare-High
            dist_int_factor(mask) = 1;
            mag_int_factor(mask) = 2;
        else  % Common-High
            dist_int_factor(mask) = 2;
            mag_int_factor(mask) = 2;
        end
        % Determine SPE factor (Expected=1, Unexpected=2)
        if mod(i_cond, 2) == 1  % Odd indices are Expected
            spe_int_factor(mask) = 1;
        else  % Even indices are Unexpected
            spe_int_factor(mask) = 2;
        end
    end

    valid_int = ~isnan(dist_int_factor) & ~isnan(spe_int_factor) & ~isnan(mag_int_factor);

    % Run 3-way ANOVA at each bin
    for i_bin = 1:n_bins
        pupil_bin = binned_pupil(:, i_bin);
        valid_idx = valid_int & ~isnan(pupil_bin);

        y = pupil_bin(valid_idx);
        g_dist = dist_int_factor(valid_idx);
        g_spe = spe_int_factor(valid_idx);
        g_mag = mag_int_factor(valid_idx);

        if length(unique(g_dist)) >= 2 && length(unique(g_spe)) >= 2 && ...
                length(unique(g_mag)) >= 2
            try
                [~, tbl, ~] = anovan(y, {g_dist, g_spe, g_mag}, ...
                    'model', 'full', ...
                    'varnames', {'Dist', 'SPE', 'Mag'}, ...
                    'display', 'off');
                % anovan output for 3-way full model:
                % Row 2: Dist main effect
                % Row 3: SPE main effect
                % Row 4: Mag main effect
                % Row 5: Dist*SPE
                % Row 6: Dist*Mag
                % Row 7: SPE*Mag
                % Row 8: Dist*SPE*Mag
                results.interaction.p_dist(i_bin) = tbl{2, 7};
                results.interaction.p_spe(i_bin) = tbl{3, 7};
                results.interaction.p_mag(i_bin) = tbl{4, 7};
                results.interaction.p_dist_spe(i_bin) = tbl{5, 7};
                results.interaction.p_dist_mag(i_bin) = tbl{6, 7};
                results.interaction.p_spe_mag(i_bin) = tbl{7, 7};
                results.interaction.p_three_way(i_bin) = tbl{8, 7};
            catch
                % Leave as NaN
            end
        end

        % Compute means and CI for each condition
        for i_cond = 1:8
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
else
    % No-AV session: set interaction results to empty/NaN
    results.interaction.trial_counts = nan(1, 8);
    results.interaction.labels = {};
    results.interaction.p_dist = nan(1, n_bins);
    results.interaction.p_spe = nan(1, n_bins);
    results.interaction.p_mag = nan(1, n_bins);
    results.interaction.p_dist_spe = nan(1, n_bins);
    results.interaction.p_dist_mag = nan(1, n_bins);
    results.interaction.p_spe_mag = nan(1, n_bins);
    results.interaction.p_three_way = nan(1, n_bins);
    results.interaction.means = nan(8, n_bins);
    results.interaction.ci_95 = nan(8, n_bins, 2);
end

%% Print Summary
fprintf('Pupil PE analysis complete.\n');
fprintf('\n--- Panel A: RPE (Magnitude-Matched) ---\n');
fprintf('  Trial counts: Rare-Low=%d, Common-Low=%d, Rare-High=%d, Common-High=%d\n', ...
    results.rpe.trial_counts(1), results.rpe.trial_counts(2), ...
    results.rpe.trial_counts(3), results.rpe.trial_counts(4));

if is_av_session
    fprintf('\n--- Panel B: SPE ---\n');
    fprintf('  Trial counts: No-Flicker=%d, Omitted=%d, Surprising=%d, ', ...
        results.spe.trial_counts(1), results.spe.trial_counts(2), ...
        results.spe.trial_counts(3));
    fprintf('Certain=%d\n', results.spe.trial_counts(4));

    fprintf('\n--- Panels C/D: RPE x SPE Interaction (8 conditions) ---\n');
    fprintf('  Trial counts:\n');
    fprintf('                     Expected  Unexpected\n');
    fprintf('    Rare-Low:        %4d      %4d\n', ...
        results.interaction.trial_counts(1), results.interaction.trial_counts(2));
    fprintf('    Common-Low:      %4d      %4d\n', ...
        results.interaction.trial_counts(3), results.interaction.trial_counts(4));
    fprintf('    Rare-High:       %4d      %4d\n', ...
        results.interaction.trial_counts(5), results.interaction.trial_counts(6));
    fprintf('    Common-High:     %4d      %4d\n', ...
        results.interaction.trial_counts(7), results.interaction.trial_counts(8));
else
    fprintf('\n--- No-AV session: SPE and Interaction analyses skipped ---\n');
end

end
