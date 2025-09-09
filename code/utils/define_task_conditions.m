function conditions = define_task_conditions(trialInfo, eventTimes, codes)
% DEFINE_TASK_CONDITIONS - Creates a struct of logical masks for trial conditions.
%
% This function first identifies valid 'tokens' task trials and then
% generates logical masks for various conditions within those trials. The
% output masks are therefore the same length as the number of tokens trials.
%
% INPUTS:
%   trialInfo  - Struct with trial-by-trial information for the whole session.
%   eventTimes - Struct with event times for each trial for the whole session.
%   codes      - Struct with unique task codes.
%
% OUTPUTS:
%   conditions - Struct with logical masks for different trial conditions,
%                filtered for tokens task trials.

% Identify valid tokens trials (task code match and reward delivered)
is_tokens_trial = (trialInfo.taskCode == codes.uniqueTaskCode_tokens) & ...
    ~cellfun(@isempty, eventTimes.rewardCell);

% Filter all relevant data structures to include only tokens trials.
% This ensures that all generated masks are of the correct length.
trialInfo.cueFile = trialInfo.cueFile(is_tokens_trial);
trialInfo.dist = trialInfo.dist(is_tokens_trial);
trialInfo.rewardAmt = trialInfo.rewardAmt(is_tokens_trial);
trialInfo.isAVTrial = trialInfo.isAVTrial(is_tokens_trial);
eventTimes.reward = eventTimes.reward(is_tokens_trial);

% A. Foundational Conditions
conditions.is_familiar = contains(trialInfo.cueFile, 'fam');
conditions.is_novel = contains(trialInfo.cueFile, 'nov');
conditions.is_normal_dist = trialInfo.dist == 1;
conditions.is_uniform_dist = trialInfo.dist == 2;
conditions.is_rewarded = eventTimes.reward > 0;

% B. Reward Magnitude / RPE Conditions
conditions.is_norm_rare_low = conditions.is_normal_dist & ...
    (trialInfo.rewardAmt <= 2);
conditions.is_norm_common = conditions.is_normal_dist & ...
    (trialInfo.rewardAmt >= 4 & trialInfo.rewardAmt <= 6);
conditions.is_norm_rare_high = conditions.is_normal_dist & ...
    (trialInfo.rewardAmt >= 8);
conditions.is_unif_low = conditions.is_uniform_dist & ...
    (trialInfo.rewardAmt <= 3);
conditions.is_unif_mid = conditions.is_uniform_dist & ...
    (trialInfo.rewardAmt >= 4 & trialInfo.rewardAmt <= 6);
conditions.is_unif_high = conditions.is_uniform_dist & ...
    (trialInfo.rewardAmt >= 7);

% C. Sensory Prediction Error (SPE) Conditions
conditions.is_flicker_certain = contains(trialInfo.cueFile, ...
    '_03.jpg') & trialInfo.isAVTrial == true;
conditions.is_flicker_surprising = contains(trialInfo.cueFile, ...
    '_02.jpg') & trialInfo.isAVTrial == true;
conditions.is_flicker_omitted = contains(trialInfo.cueFile, ...
    '_02.jpg') & trialInfo.isAVTrial == false;
conditions.is_noflicker_certain = contains(trialInfo.cueFile, ...
    '_01.jpg');

% D. Key RPE x SPE Interaction Conditions
conditions.is_common_reward_no_spe = conditions.is_norm_common & ...
    conditions.is_noflicker_certain;
conditions.is_rare_high_reward_no_spe = conditions.is_norm_rare_high & ...
    conditions.is_noflicker_certain;
conditions.is_common_reward_with_spe = conditions.is_norm_common & ...
    conditions.is_flicker_surprising;
conditions.is_rare_high_reward_with_spe = conditions.is_norm_rare_high & ...
    conditions.is_flicker_surprising;

end
