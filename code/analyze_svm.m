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
%   core_data:  Struct containing the binned firing rates and time vectors.
%               - .binned_spikes: [n_trials x n_neurons x n_time_bins]
%               - .time_vector:   [1 x n_time_bins]
%   conditions: Struct of logical masks for different trial conditions.
%   svm_plan:   A single element from the svm_plan struct, defining the
%               specific comparison to be made.
%
% OUTPUTS:
%   svm_results: Struct containing the analysis results.
%                - .accuracy:     [1 x n_time_bins] classification accuracy
%                - .accuracy_ci:  [2 x n_time_bins] 95% CI for accuracy
%                - .time_vector:  [1 x n_time_bins] copied from core_data

%% Setup Paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% --- 1. Prepare Data ---
% Get the condition masks from the plan
cond1_mask = conditions.(svm_plan.cond1);
cond2_mask = conditions.(svm_plan.cond2);

% Combine masks to select all relevant trials
combined_mask = cond1_mask | cond2_mask;

% Filter the spike data to include only the trials for this comparison
binned_spikes = core_data.binned_spikes(combined_mask, :, :);

% Create the label vector Y
% '1' for cond1, '2' for cond2
Y = cond2_mask(combined_mask) + 1;

% Get dimensions
[~, n_neurons, n_time_bins] = size(binned_spikes);

% Pre-allocate results structure
svm_results.accuracy = nan(1, n_time_bins);
svm_results.accuracy_ci = nan(2, n_time_bins);
svm_results.time_vector = core_data.time_vector;

%% --- 2. Time-Resolved SVM Analysis ---
% Loop over each time bin to perform the classification
for t = 1:n_time_bins
    % Prepare the data matrix X for the current time bin
    % X should be [n_trials x n_neurons]
    X = squeeze(binned_spikes(:, :, t));

    % If there's no variance in the data for this bin, skip it
    if range(X(:)) == 0
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
