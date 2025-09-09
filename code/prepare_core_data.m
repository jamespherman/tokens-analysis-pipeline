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

%% Define Alignment Events and Trial Indices
% Define the list of events to align the data to
alignment_events = {'CUE_ON', 'outcomeOn', 'reward'};

% Load the structure containing all task codes
codes = initCodes();

% Find the indices of trials that belong to the 'tokens' task
tokens_trial_indices = find(session_data.trialInfo.taskCode == ...
    codes.uniqueTaskCode_tokens);

%% Prepare Neuronal Data
% Pass the alignment events to the neuronal data preparation function
core_data.spikes = prepare_neuronal_data(session_data, ...
    selected_neurons, tokens_trial_indices, alignment_events);

%% Prepare Pupil Data
% Pass the alignment events to the pupil data preparation function
core_data.pupil = prepare_pupil_data(session_data, ...
    tokens_trial_indices, alignment_events);

end
