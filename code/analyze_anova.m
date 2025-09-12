%% analyze_anova.m
%
% Performs an N-way ANOVA on firing rate data for each neuron and time bin
% across three alignment events. It tests for the main effects and pairwise
% interactions of RPE, flicker condition, and reward distribution.
%
% INPUTS:
%   session_data - The main data structure for the session.
%   core_data    - A struct containing binned firing rates.
%   conditions   - A struct with logical masks for trial conditions.
%
% OUTPUT:
%   session_data - The input session_data struct with a new field,
%                  'analysis.anova_results', containing the p-values.
%
% Author: Jules
% Date: 2025-09-11
%

function session_data = analyze_anova(session_data, core_data, conditions)

%% Setup Paths
% Add the 'utils' directory to the path so that helper functions can be
% found.
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Define Alignment Events
alignment_events = {'CUE_ON', 'outcomeOn', 'reward'};

%% Prepare Factors for ANOVA
% These factors are constant across all neurons and time bins. They are
% prepared once here to avoid redundant calculations inside the main loops.

% To get the correct trial-by-trial data, we must first identify the same
% set of rewarded 'tokens' trials that were used to create the core_data
% and conditions structs.
codes = initCodes();
tempCueFile = session_data.trialInfo.cueFile;
tempCueFile(cellfun(@isempty, tempCueFile)) = {''};

% define one variable indexing all tokens task trials including those that
% are not going to be part of this anova analysis:
is_tokens_trial_all = (session_data.trialInfo.taskCode == ...
    codes.uniqueTaskCode_tokens) & ...
    ~cellfun(@isempty, session_data.eventTimes.rewardCell);

% define another variable that indexes only the tokens task trials that
% will be used for the anova analysis (those that had a cue image):
is_tokens_trial = is_tokens_trial_all & ...
    ~contains(tempCueFile, 'blank');
n_trials = sum(is_tokens_trial);

% Factor 1: RPE (Continuous Predictor)
% Calculated as the reward amount on each trial minus the mean reward amount
% across all valid tokens trials.
rewardAmt = session_data.trialInfo.rewardAmt(is_tokens_trial);
rpe_predictor = rewardAmt - mean(rewardAmt);

% Factor 2: Distribution (Categorical, 2 Levels: 'Normal', 'Uniform')
% A cell array of strings defining the reward distribution for each trial.
% first make a logical index of the same size as
% 'conditions.is_normal_dist' / 'conditions.is_uniform_dist' so we can use
% those to define 'dist_factor':
gDist = find(~contains(tempCueFile(is_tokens_trial_all), 'blank'));
dist_factor = cell(n_trials, 1);
dist_factor(conditions.is_normal_dist(gDist)) = {'Normal'};
dist_factor(conditions.is_uniform_dist(gDist)) = {'Uniform'};

% Factor 3: Flicker (Categorical, 4 Levels)
% A cell array of strings defining the flicker condition for each trial.
flicker_factor = cell(n_trials, 1);
flicker_factor(conditions.is_noflicker_certain(gDist)) = {'CertainAbsent'};
flicker_factor(conditions.is_flicker_certain(gDist)) = {'CertainPresent'};
flicker_factor(conditions.is_flicker_omitted(gDist)) = {'UncertainAbsent'};
flicker_factor(conditions.is_flicker_surprising(gDist)) = {'UncertainPresent'};


%% Main Analysis Loop
for i_event = 1:numel(alignment_events)
    event_name = alignment_events{i_event};

    % Get dimensions from core_data for this event
    [n_neurons, n_trials, n_bins] = size(...
        core_data.spikes.(event_name).rates);

    % Pre-allocate the results structure for the current event. This improves
    % performance and code clarity. The p-values for each model term will be
    % stored in a matrix of size [n_neurons x n_bins].
    term_names = {'rpe', 'dist', 'flicker', 'rpe_dist', 'rpe_flicker', ...
        'dist_flicker'};
    for i_term = 1:numel(term_names)
        p_field_name = ['p_' term_names{i_term}];
        results.(event_name).(p_field_name) = nan(n_neurons, n_bins);
    end

    % Loop through each neuron and each time bin to perform the ANOVA
    for i_neuron = 1:n_neurons
        for i_bin = 1:n_bins

            % Extract the firing rates for the current neuron and time bin
            % across all trials. Squeeze removes singleton dimensions.
            firing_rates = squeeze(core_data.spikes.(event_name).rates( ...
                i_neuron, gDist, i_bin));

            % Anovan requires at least two groups for each categorical
            % variable. If NaNs in the firing rate data cause a factor to
            % have only one level, anovan will error. We check for this
            % and skip the bin if the condition isn't met.
            valid_trials = ~isnan(firing_rates);
            if numel(unique(dist_factor(valid_trials))) < 2 || ...
               numel(unique(flicker_factor(valid_trials))) < 2
                continue;
            end

            % Perform the N-way ANOVA. Use a try-catch block to gracefully
            % handle cases where anovan might fail (e.g., rank deficiency).
            try
                [~, tbl, ~] = anovan(firing_rates, ...
                    {rpe_predictor, dist_factor, flicker_factor}, ...
                    'model', 'interaction', ...
                    'continuous', 1, ...
                    'varnames', {'rpe', 'dist', 'flicker'}, ...
                    'display', 'off');

                % p-values are in the 7th column of the output table.
                % We store the p-value for each main effect and interaction.
                results.(event_name).p_rpe(i_neuron, i_bin) = tbl{2, 7};
                results.(event_name).p_dist(i_neuron, i_bin) = tbl{3, 7};
                results.(event_name).p_flicker(i_neuron, i_bin) = tbl{4, 7};
                results.(event_name).p_rpe_dist(i_neuron, i_bin) = tbl{5, 7};
                results.(event_name).p_rpe_flicker(i_neuron, i_bin) = tbl{6, 7};
                results.(event_name).p_dist_flicker(i_neuron, i_bin) = tbl{7, 7};

            catch ME
                % If anovan fails, the p-values for this bin remain NaN.
                % This is expected for bins with insufficient data.
                % Uncomment the line below for detailed debugging.
                % fprintf('Could not run ANOVA for neuron %d, bin %d, event %s: %s\n', ...
                %     i_neuron, i_bin, event_name, ME.message);
            end

        end % End of time bin loop
    end % End of neuron loop
end % End of event loop

%% Finalize Output
% Assign the completed results structure to a new field in the session_data
% struct. This modified struct is the final output of the function.
session_data.analysis.anova_results = results;

fprintf('N-way ANOVA analysis complete.\n');

end
