%% prepare_pupil_data.m
%
%   Generates preprocessed, event-aligned pupil data traces.
%
% INPUTS:
%   session_data         - The main data structure.
%   tokens_trial_indices - Indices for trials of the 'tokens' task.
%
% OUTPUT:
%   aligned_pupil        - A struct with aligned pupil data.
%
% Author: Jules
% Date: 2025-09-08

function aligned_pupil = prepare_pupil_data(session_data, ...
    tokens_trial_indices)

%% Define Alignment and Preprocessing Parameters
alignment_events = {'cueOn', 'outcomeOn'};
time_window = [-0.5, 1.5]; % in seconds
sample_rate = session_data.display.fr; % Hz, from display refresh rate
n_samples = round((time_window(2) - time_window(1)) * sample_rate);
time_vector = linspace(time_window(1), time_window(2), n_samples);

baseline_window = [-0.5, -0.1]; % in seconds, relative to event
deblink_threshold = 0.1; % Pupil size values below this are artifacts
smoothing_window_ms = 50; % ms
smoothing_window_samples = round(smoothing_window_ms / 1000 * sample_rate);

n_tokens_trials = numel(tokens_trial_indices);

%% Process Each Alignment Event
for i_event = 1:numel(alignment_events)
    event_name = alignment_events{i_event};

    % Initialize storage for this event
    aligned_traces = nan(n_tokens_trials, n_samples);

    % Get the event times for the current alignment event
    event_times = session_data.behavior.eventTimes.(event_name)...
        (tokens_trial_indices);

    % Loop through each trial
    for i_trial = 1:n_tokens_trials
        trial_idx = tokens_trial_indices(i_trial);

        % a. Check for Data
        if isempty(session_data.pupil.raw{trial_idx}) || ...
           isempty(session_data.pupil.t{trial_idx})
            continue;
        end

        pupil_trace = session_data.pupil.raw{trial_idx};
        pupil_time = session_data.pupil.t{trial_idx} - event_times(i_trial);

        % b. De-blink
        pupil_trace(pupil_trace < deblink_threshold) = nan;

        % c. Smooth
        smoothed_trace = movmean(pupil_trace, smoothing_window_samples, ...
            'omitnan');

        % d. Baseline Correction
        baseline_indices = pupil_time >= baseline_window(1) & ...
                           pupil_time <= baseline_window(2);
        baseline_mean = mean(smoothed_trace(baseline_indices), 'omitnan');

        if isnan(baseline_mean) || baseline_mean == 0
            continue; % Cannot normalize if baseline is invalid
        end

        % e. Normalization
        normalized_trace = (smoothed_trace - baseline_mean) / baseline_mean;

        % f. Align & Store
        % Interpolate to a common time base
        aligned_traces(i_trial, :) = interp1(pupil_time, ...
            normalized_trace, time_vector, 'linear', nan);
    end

    % Store the final matrix and its time vector
    aligned_pupil.(event_name).traces = aligned_traces;
    aligned_pupil.(event_name).time_vector = time_vector;
end

end
