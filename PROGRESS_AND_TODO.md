# Tokens AV: Immediate Analysis Priorities (Pre-ACNP)

## Goal
Generate 1-2 figures showing behavioral/pupil evidence that primates show AV-modulated responses in the tokens task, to give Vinny something concrete to show Holly at ACNP (leaving Sunday, ~1 week).

---

## Priority 1: Pupil Time-Varying Analysis

### Data Requirements
- Pupil traces from tokens_AV sessions
- Trial event times (cue onset, outcome onset, reward onset)
- Condition labels (Distribution × AV Probability)
- Token counts per trial

### Preprocessing Steps
```matlab
% 1. Load pupil data
% 2. Detect blinks (pupil = 0 or rapid changes)
% 3. Interpolate blinks (linear or cubic spline)
% 4. Baseline normalize: subtract mean of -200 to 0ms pre-cue
% 5. Z-score within session (optional, helps cross-session comparison)
% 6. Reject trials with >30% missing data
```

### Analysis 1A: Condition Averages
```matlab
% Align to outcome onset
% Average pupil trace by condition (2 dist × 3 AV prob = 6 conditions)
% Plot with SEM shading
% Time window: -500ms to +2000ms relative to outcome
```

**Expected figure**: 6 traces showing whether AV conditions have larger/different pupil responses

### Analysis 1B: Time-Varying ANOVA
```matlab
% At each 100ms bin from outcome onset to +1500ms:
for t = 1:15  % 100ms bins
    % Extract pupil value in this bin for each trial
    pupil_bin = mean(pupil_aligned(:, bin_start:bin_end), 2);
    
    % 2-way ANOVA: Distribution × AV_Probability
    [p_dist(t), p_av(t), p_interaction(t)] = anovan(pupil_bin, ...
        {distribution, av_probability}, 'display', 'off');
end

% Plot: 3 lines showing p-values over time
% Highlight significant epochs (p < 0.05, corrected)
```

**Key question**: Is there an INTERACTION between AV probability and distribution? This would suggest AV isn't just an additive arousal effect.

### Analysis 1C: AV Surprise (50% condition only)
```matlab
% For 50% AV cues only:
% Compare trials where AV was presented vs. not presented
idx_50 = av_probability == 0.5;

% t-test or mixed model at each time bin
for t = 1:15
    pupil_av_present = pupil_bin(idx_50 & av_present);
    pupil_av_absent = pupil_bin(idx_50 & ~av_present);
    [~, p_av_surprise(t)] = ttest2(pupil_av_present, pupil_av_absent);
end
```

**Expected result**: Larger pupil response when AV is presented vs. absent (within the 50% condition)

---

## Priority 2: Reward Surprise × AV Interaction

### Compute Reward Surprise
```matlab
% For Normal distribution (mean = 5):
reward_surprise_norm = token_count - 5;

% For Uniform distribution (mean = 5.5):
reward_surprise_uni = token_count - 5.5;

% Combine into single variable
reward_surprise = nan(n_trials, 1);
reward_surprise(distribution == 'Normal') = reward_surprise_norm;
reward_surprise(distribution == 'Uniform') = reward_surprise_uni;
```

### Analysis 2A: Pupil × Reward Surprise Regression
```matlab
% For a specific time window (e.g., 300-800ms post-outcome):
pupil_response = mean(pupil_aligned(:, 300:800ms), 2);

% Linear model
mdl = fitlme(tbl, 'pupil_response ~ reward_surprise * av_present + (1|session)');

% Key terms:
% - reward_surprise: Does pupil track reward magnitude?
% - av_present: Does AV elevate pupil overall?
% - reward_surprise:av_present: INTERACTION - does AV change reward sensitivity?
```

### Analysis 2B: Visualization
```matlab
% Scatter plot: reward_surprise (x) vs. pupil_response (y)
% Color code by AV present (yes/no)
% Add regression lines for each

% OR: Bin reward surprise into quintiles, plot pupil by bin, separate lines for AV/no-AV
```

---

## Priority 3: Task Schematic Figure

### Elements Needed
1. Trial timeline with event markers
2. Example cue images (or schematic representations)
3. Token display (with and without AV flicker)
4. Condition matrix showing 2×3 design

### Timeline Schematic
```
ITI → Cue → Fixation Hold → Outcome (±AV) → Cash-in → End
       ↓         ↓              ↓
    500ms     500ms          1000ms → reward pulses
```

---

## Deliverable Checklist

### Minimum for ACNP
- [ ] Task schematic (can be hand-drawn or PowerPoint)
- [ ] Pupil traces by AV probability (showing main effect of AV)
- [ ] One statistical test showing interaction or AV surprise effect

### Ideal for ACNP
- [ ] Full 2×3 pupil analysis with ANOVA statistics
- [ ] Reward surprise × AV interaction plot
- [ ] Preliminary SC sensory PE result (if time permits)

---

## Quick Sanity Checks

Before diving deep, verify:

1. **Do you have usable pupil data?**
   - Check a few sessions for data quality
   - Confirm alignment to events is correct

2. **How many trials per condition?**
   - Need ~50+ per condition for stable estimates
   - Check if pooling across sessions is needed

3. **Which monkeys/sessions?**
   - Document which animals contribute to analysis
   - Note any exclusions

---

## Code Snippets

### Loading and Aligning Pupil Data (pseudocode)
```matlab
function [pupil_aligned, trial_info] = load_pupil_session(session_path)
    % Load raw pupil
    load(fullfile(session_path, 'pupil_data.mat'));
    
    % Load trial events
    load(fullfile(session_path, 'trial_events.mat'));
    
    % Define alignment windows
    pre_outcome = 500;  % ms before outcome
    post_outcome = 2000;  % ms after outcome
    
    % Loop through successful trials
    for tr = 1:n_trials
        if trial_info(tr).outcome == 'success'
            outcome_time = trial_events(tr).outcome_onset;
            
            % Extract pupil around outcome
            t_start = outcome_time - pre_outcome;
            t_end = outcome_time + post_outcome;
            
            pupil_aligned(tr, :) = interp_pupil(pupil_raw, t_start:t_end);
        end
    end
end
```

### Quick ANOVA Wrapper
```matlab
function [p_main1, p_main2, p_interaction] = quick_anova_2way(y, factor1, factor2)
    % Wrapper for 2-way ANOVA returning just the p-values
    [p, tbl] = anovan(y, {factor1, factor2}, ...
        'model', 'interaction', ...
        'varnames', {'Factor1', 'Factor2'}, ...
        'display', 'off');
    
    p_main1 = p(1);
    p_main2 = p(2);
    p_interaction = p(3);
end
```

---

## Notes for James

1. **Start with what you have**: Even partial results are useful. If pupil data is messy, show what you can.

2. **Focus on the 50% AV condition**: This is where the magic is - it's the only condition with genuine AV surprise.

3. **Simple is fine**: A clean 2×3 mean plot with error bars and one significant statistic is better than complex analyses with unclear results.

4. **Backup plan**: If pupil doesn't work out, can you show the SC sensory PE finding (unexpected > expected AV) as a teaser?

---

*Created: January 5, 2026*