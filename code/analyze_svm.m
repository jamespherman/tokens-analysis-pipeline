%% analyze_svm.m
%
% Performs a time-resolved, cross-validated SVM classification analysis.
%
% This function takes core data, trial conditions, and a specific
% comparison plan, then runs a linear SVM classifier for each time bin. It
% calculates classification accuracy and 95% confidence intervals.
%
% Author: Jules
% Date: 2025-09-15
%

function svm_results = analyze_svm(core_data, conditions, svm_plan)
% ANALYZE_SVM Time-resolved SVM classification
%
% INPUTS:
%   core_data:  Struct containing processed data. It is expected to have
%               the following substructure:
%               - .spikes.(eventName).rates: [n_neurons x n_trials x n_time_bins]
%               - .spikes.(eventName).time_vector: [1 x n_time_bins]
%   conditions: Struct of logical masks for different trial conditions.
%   svm_plan:   A single element from the svm_plan struct, defining the
%               specific comparison to be made. Must contain an `.event`
%               field specifying which eventName to use from core_data.
%
% OUTPUTS:
%   svm_results: Struct containing the analysis results.
%                - .accuracy:     [1 x n_time_bins] classification accuracy
%                - .accuracy_ci:  [2 x n_time_bins] 95% CI for accuracy
%                - .time_vector:  [1 x n_time_bins] copied from core_data's
%                                 event-specific time vector.

%% Setup Paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% --- 1. Prepare Data ---
% Get the condition masks from the plan
cond1_mask = conditions.(svm_plan.cond1);
cond2_mask = conditions.(svm_plan.cond2);

% Combine masks to select all relevant trials
combined_mask = cond1_mask | cond2_mask;

% Get the correct binned firing rates using the event from the plan
% The raw data is [n_neurons x n_trials x n_time_bins], so we need to
% permute it to [n_trials x n_neurons x n_time_bins] for the analysis.
event_rates = core_data.spikes.(svm_plan.event).rates;
event_rates_permuted = permute(event_rates, [2, 1, 3]);

% Filter the spike data to include only the trials for this comparison
binned_spikes = event_rates_permuted(combined_mask, :, :);

% Create the label vector Y
% '1' for cond1, '2' for cond2
Y = cond2_mask(combined_mask) + 1;

% Get dimensions
[~, n_neurons, n_time_bins] = size(binned_spikes);

% Pre-allocate results structure
svm_results.accuracy = nan(1, n_time_bins);
svm_results.accuracy_ci = nan(2, n_time_bins);
% Get the time vector for the specific event alignment
svm_results.time_vector = core_data.spikes.(svm_plan.event).time_vector;

%% --- 2. Time-Resolved SVM Analysis ---
% Loop over each time bin to perform the classification
for t = 1:n_time_bins

    % Create the label vector Y
    % '1' for cond1, '2' for cond2
    Y = cond2_mask(combined_mask) + 1;

    % Prepare the data matrix X for the current time bin
    % X should be [n_trials x n_neurons]
    X = squeeze(binned_spikes(:, :, t));

    % get rid of rows that are all NaN:
    allNanRows = all(isnan(X), 2);
    X(allNanRows,:) = [];
    Y(allNanRows)   = [];

    % If there's no variance in the data for this bin, skip it
    xRange = range(X(:));
    if xRange == 0 || isempty(xRange)
        continue;
    end

    % Train a linear SVM with 10-fold cross-validation
    % Using 'Prior', 'Uniform' to handle unbalanced classes
    SVMModel = fitcsvm(X, Y, 'KernelFunction', 'linear', ...
        'Standardize', true, 'CrossVal', 'on', 'Prior', 'Uniform');

    % Get the cross-validated predictions
    predictions = kfoldPredict(SVMModel);

    % Calculate accuracy
    n_correct = sum(predictions == Y);
    n_total = numel(Y);

    % Use binofit to get accuracy and 95% confidence intervals
    [acc, acc_ci] = binofit(n_correct, n_total);

    % Store the results
    svm_results.accuracy(t) = acc;
    svm_results.accuracy_ci(:, t) = acc_ci;
end

end
