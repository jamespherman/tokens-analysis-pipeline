function analysis_plan = define_analysis_plan()
% DEFINE_ANALYSIS_PLAN Central configuration for tokens analysis pipeline
%
%   analysis_plan = DEFINE_ANALYSIS_PLAN()
%   Creates and returns a struct that defines all analyses and comparisons
%   to be run in the tokens analysis pipeline. This function programmatically
%   generates the analysis plan by calling define_task_conditions to get the
%   canonical condition names.
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
%   See also ANALYZE_BASELINE_COMPARISON, ANALYZE_ROC_COMPARISON, DEFINE_TASK_CONDITIONS

% Call define_task_conditions to get the canonical condition definitions.
% We call it without inputs, as we only need the condition_defs struct.
[~, ~, condition_defs] = define_task_conditions();

% Programmatically build the "Baseline vs. Post-Event Activity" analysis plan
analysis_plan.baseline_comparison.conditions_to_run = [ ...
    condition_defs.RPE_normal, ...
    condition_defs.SPE ...
    ];

% Programmatically build the "Bin-by-Bin ROC Comparison" analysis plan
comparisons_def = { ...
    'Dist_at_Cue',    'CUE_ON',    condition_defs.distribution; ...
    'RPE_at_Outcome', 'outcomeOn', condition_defs.RPE_comparison_pair; ...
    'RPE_at_Reward',  'reward',    condition_defs.RPE_comparison_pair; ...
    'SPE_at_Outcome', 'outcomeOn', condition_defs.SPE_comparison_pair ...
};

comparisons = struct('name', {}, 'event', {}, 'cond1', {}, 'cond2', {});
for i = 1:size(comparisons_def, 1)
    comparisons(i).name  = comparisons_def{i, 1};
    comparisons(i).event = comparisons_def{i, 2};
    comparisons(i).cond1 = comparisons_def{i, 3}{1};
    comparisons(i).cond2 = comparisons_def{i, 3}{2};
end

analysis_plan.roc_comparison.comparisons_to_run = comparisons;


% "N-way ANOVA" Analysis
analysis_plan.anova.run = true;
analysis_plan.anova.fields_to_aggregate = { ...
    'p_value_reward', 'p_value_stim_id', 'p_value_interaction', ...
    'p_value_flicker', 'p_value_flicker_x_reward' ...
    };

end
