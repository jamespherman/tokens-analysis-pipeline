function analysis_plan = define_analysis_plan()
% DEFINE_ANALYSIS_PLAN Central configuration for tokens analysis pipeline
%
%   analysis_plan = DEFINE_ANALYSIS_PLAN()
%   Creates and returns a struct that defines all analyses and comparisons
%   to be run in the tokens analysis pipeline.
%
%   OUTPUT
%   - analysis_plan: A struct with the following fields:
%       .baseline_comparison
%           .conditions_to_run: Cell array of condition names for
%                               "Baseline vs. Post-Event Activity" analysis.
%       .roc_comparison
%           .comparisons_to_run: Structure array defining bin-by-bin ROC
%                                comparisons. Each element has the fields:
%                                .name, .event, .cond1, .cond2
%
%   See also ANALYZE_BASELINE_COMPARISON, ANALYZE_ROC_COMPARISON

% "Baseline vs. Post-Event Activity" Analysis
analysis_plan.baseline_comparison.conditions_to_run = { ...
    'is_normal_dist', ...
    'is_uniform_dist', ...
    'is_norm_common', ...
    'is_norm_rare_high', ...
    'is_flicker_surprising', ...
    'is_flicker_certain' ...
    };

% "Bin-by-Bin ROC Comparison" Analysis
analysis_plan.roc_comparison.comparisons_to_run = struct( ...
    'name',  { 'Dist_at_Cue',        'RPE_at_Outcome',      'RPE_at_Reward',      'SPE_at_Outcome' }, ...
    'event', { 'CUE_ON',             'outcomeOn',           'reward',             'outcomeOn' }, ...
    'cond1', { 'is_normal_dist',     'is_norm_common',      'is_norm_common',     'is_flicker_certain' }, ...
    'cond2', { 'is_uniform_dist',    'is_norm_rare_high',   'is_norm_rare_high',  'is_flicker_surprising' } ...
    );

% "N-way ANOVA" Analysis
analysis_plan.anova.run = true;

end
