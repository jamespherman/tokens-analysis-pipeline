function [conditions, is_av_session] = define_task_conditions(trialInfo, eventTimes, unique_id)

% DEFINE_TASK_CONDITIONS - Creates a struct of logical masks for trial conditions.
%
% This function first identifies valid 'tokens' task trials and then
% generates logical masks for various conditions within those trials. The
% output masks are therefore the same length as the number of tokens trials.
%
% INPUTS:
%   trialInfo  - Struct with trial-by-trial information for the whole session.
%   eventTimes - Struct with event times for each trial for the whole session.
%   unique_id  - String, a unique identifier for the session, used for saving figures.
%
% OUTPUTS:
%   conditions - Struct with logical masks for different trial conditions,
%                filtered for tokens task trials.
%   is_av_session - Logical flag indicating if the session includes AV trials.

codes = initCodes;

% Determine if this is an AV session by checking for the 'isAVTrial' field.
is_av_session = isfield(trialInfo, 'isAVTrial');

% Identify valid tokens trials (task code match and reward delivered)
is_tokens_trial = (trialInfo.taskCode == codes.uniqueTaskCode_tokens) & ...
    ~cellfun(@isempty, eventTimes.rewardCell) & ...
    ~cellfun(@isempty, trialInfo.cueFile);

% Filter all relevant data structures to include only tokens trials.
% This ensures that all generated masks are of the correct length.
trialInfo.cueFile = trialInfo.cueFile(is_tokens_trial);
trialInfo.dist = trialInfo.dist(is_tokens_trial);
trialInfo.rewardAmt = trialInfo.rewardAmt(is_tokens_trial);
if is_av_session
    trialInfo.isAVTrial = trialInfo.isAVTrial(is_tokens_trial);
end
eventTimes.reward = eventTimes.reward(is_tokens_trial);

% A. Foundational Conditions
conditions.is_familiar = contains(trialInfo.cueFile, 'fam');
conditions.is_novel = contains(trialInfo.cueFile, 'nov');
conditions.is_normal_dist = trialInfo.dist == 1;
conditions.is_uniform_dist = trialInfo.dist == 2;
conditions.is_rewarded = eventTimes.reward > 0;

% B. Reward Magnitude / RPE Conditions
% For the normal distribution, calculate dynamic thresholds
reward_norm = trialInfo.rewardAmt(conditions.is_normal_dist);
norm_thresholds = prctile(reward_norm, [25 75]);
norm_p25 = norm_thresholds(1);
norm_p75 = norm_thresholds(2);

conditions.is_norm_rare_low = conditions.is_normal_dist & ...
    (trialInfo.rewardAmt <= norm_p25);
conditions.is_norm_common = conditions.is_normal_dist & ...
    (trialInfo.rewardAmt > norm_p25 & trialInfo.rewardAmt < norm_p75);
conditions.is_norm_rare_high = conditions.is_normal_dist & ...
    (trialInfo.rewardAmt >= norm_p75);

% For the uniform distribution, calculate dynamic thresholds
reward_unif = trialInfo.rewardAmt(conditions.is_uniform_dist);
unif_thresholds = prctile(reward_unif, [25 75]);
unif_p25 = unif_thresholds(1);
unif_p75 = unif_thresholds(2);

conditions.is_unif_low = conditions.is_uniform_dist & ...
    (trialInfo.rewardAmt <= unif_p25);
conditions.is_unif_mid = conditions.is_uniform_dist & ...
    (trialInfo.rewardAmt > unif_p25 & trialInfo.rewardAmt < unif_p75);
conditions.is_unif_high = conditions.is_uniform_dist & ...
    (trialInfo.rewardAmt >= unif_p75);

% C. Sensory Prediction Error (SPE) Conditions
if is_av_session
    conditions.is_flicker_certain = contains(trialInfo.cueFile, ...
        '_03.jpg') & trialInfo.isAVTrial == true;
    conditions.is_flicker_surprising = contains(trialInfo.cueFile, ...
        '_02.jpg') & trialInfo.isAVTrial == true;
    conditions.is_flicker_omitted = contains(trialInfo.cueFile, ...
        '_02.jpg') & trialInfo.isAVTrial == false;
    conditions.is_noflicker_certain = contains(trialInfo.cueFile, ...
        '_01.jpg');
else
    % If it's not an AV session, create all-false placeholders for SPE conditions
    % to ensure the 'conditions' struct has a consistent field structure.
    num_trials = length(conditions.is_familiar); % Get length from an existing field
    conditions.is_flicker_certain = false(num_trials, 1);
    conditions.is_flicker_surprising = false(num_trials, 1);
    conditions.is_flicker_omitted = false(num_trials, 1);
    conditions.is_noflicker_certain = false(num_trials, 1);
end

% D. Key RPE x SPE Interaction Conditions
conditions.is_common_reward_no_spe = conditions.is_norm_common & ...
    conditions.is_noflicker_certain;
conditions.is_rare_high_reward_no_spe = conditions.is_norm_rare_high & ...
    conditions.is_noflicker_certain;
conditions.is_common_reward_with_spe = conditions.is_norm_common & ...
    conditions.is_flicker_surprising;
conditions.is_rare_high_reward_with_spe = conditions.is_norm_rare_high & ...
    conditions.is_flicker_surprising;


% E. Generate and Save Diagnostic Plot
% Create a plot to visualize the reward distributions and thresholds.
fig = figure('Visible', 'off', 'Position', [100, 100, 800, 400]);

% Panel 1: Normal Distribution
mySubPlot([1, 2, 1]);
histogram(reward_norm);
title('Normal Distribution Outcomes');
xlabel('Reward Amount');
ylabel('Count');
hold on;
xline(norm_p25, '--r', 'LineWidth', 2);
xline(norm_p75, '--r', 'LineWidth', 2);
hold off;

% Panel 2: Uniform Distribution
mySubPlot([1, 2, 2]);
histogram(reward_unif);
title('Uniform Distribution Outcomes');
xlabel('Reward Amount');
hold on;
xline(unif_p25, '--r', 'LineWidth', 2);
xline(unif_p75, '--r', 'LineWidth', 2);
hold off;

% Add a main title for the whole figure
sgtitle(sprintf('Reward Distributions for Session: %s', unique_id), 'Interpreter', 'none');

% Save the figure to the 'figures' directory. The script is run from the
% 'code/' directory, so we use a relative path.
figures_dir = '../figures';
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
file_name = fullfile(figures_dir, sprintf('%s_reward_distributions.pdf', unique_id));
pdfSave(file_name, [11 8.5], fig);

% Close the figure to free up memory
close(fig);

end
