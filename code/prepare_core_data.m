%% prepare_core_data.m
%
%   This function serves as the main wrapper for preparing analysis-ready
%   data structures for the tokens task. It orchestrates the processing of
%   both neuronal and pupil data by calling specialized sub-functions.
%
% INPUTS:
%   session_data     - The main data structure containing all session info.
%   selected_neurons - A logical vector indicating which neurons to include.
%
% OUTPUT:
%   core_data        - A struct containing the processed neuronal and pupil
%                      data.
%
% Author: Jules
% Date: 2025-09-08

function core_data = prepare_core_data(session_data, selected_neurons)

%% Setup Paths
% Add the 'utils' directory to the path so that helper functions can be
% found.
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Identify Tokens Task Trials
% Load the structure containing all task codes
codes = initCodes();

% Find the indices of trials that belong to the 'tokens' task
tokens_trial_indices = find(session_data.trialInfo.taskCode == ...
    codes.uniqueTaskCode_tokens);

%% Prepare Neuronal Data
% Call the sub-function to get aligned spike rate matrices
core_data.spikes = prepare_neuronal_data(session_data, ...
    selected_neurons, tokens_trial_indices);

%% Prepare Pupil Data
% Call the sub-function to get preprocessed and aligned pupil traces
core_data.pupil = prepare_pupil_data(session_data, ...
    tokens_trial_indices);

end
