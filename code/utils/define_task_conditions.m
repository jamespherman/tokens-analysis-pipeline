function conditions = define_task_conditions(trialInfo, eventTimes)
% DEFINE_TASK_CONDITIONS - Creates a struct of logical masks for trial conditions.
%
% INPUTS:
%   trialInfo  - Struct with trial-by-trial information.
%   eventTimes - Struct with event times for each trial.
%
% OUTPUTS:
%   conditions - Struct with logical masks for different trial conditions.

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
