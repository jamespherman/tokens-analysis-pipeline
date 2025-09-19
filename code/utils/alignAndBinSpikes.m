function [binned_counts, bin_centers, aligned_spikes] = alignAndBinSpikes(varargin)
%
% [binned_counts, bin_centers, aligned_spikes] = alignAndBinSpikes(...
%   spikeTimes, eventTimes, minTime, maxTime, binWidth, 'step_size', step_size)
%
% alignAndBinSpikes: Aligns spike times to given event times and bins them.
%
% INPUTS:
%   spikeTimes:      Vector of spike times for a single unit.
%   eventTimes:      Times at which the events occur.
%   minTime, maxTime: Define the time window for alignment (seconds).
%   binWidth:        Width of the bins for the histogram (seconds).
%
% OPTIONAL INPUTS:
%   step_size:       The step size for the sliding window (seconds).
%                    If not provided, defaults to binWidth for non-overlapping bins.
%
% OUTPUTS:
%   binned_counts:   Matrix of spike counts (nEventTimes x nBins).
%   bin_centers:     Vector of bin center times (1 x nBins).
%   aligned_spikes:  Vector of all spike times aligned to their respective events.
%

% --- Input parsing ---
p = inputParser;
addRequired(p, 'spikeTimes', @isnumeric);
addRequired(p, 'eventTimes', @isnumeric);
addRequired(p, 'minTime', @isnumeric);
addRequired(p, 'maxTime', @isnumeric);
addRequired(p, 'binWidth', @isnumeric);
addOptional(p, 'step_size', -1, @isnumeric); % Use a sentinel value for default

parse(p, varargin{:});

spikeTimes = p.Results.spikeTimes;
eventTimes = p.Results.eventTimes;
minTime = p.Results.minTime;
maxTime = p.Results.maxTime;
binWidth = p.Results.binWidth;
step_size = p.Results.step_size;

% If step_size was not provided, default to binWidth
if step_size == -1
    step_size = binWidth;
end

% --- Bin Definition ---
bin_starts  = minTime:step_size:(maxTime - binWidth);
bin_ends    = bin_starts + binWidth;
bin_centers  = mean([bin_starts; bin_ends]);
nBins = length(bin_starts);
nEventTimes = length(eventTimes);

% --- Unified Spike Counting ---
binned_counts = zeros(nEventTimes, nBins);
aligned_spikes_cell = cell(nEventTimes, 1);

for i = 1:nEventTimes
    % Align spikes to the current event time
    event_window = spikeTimes > eventTimes(i) + minTime & ...
        spikeTimes < eventTimes(i) + maxTime;

    % store aligned spikes
    aligned_spikes_cell{i} = spikeTimes(event_window) - eventTimes(i);
    
    % count and store binned counts
    for j = 1:nBins
        binned_counts(i, j) = ...
            nnz(aligned_spikes_cell{i} > bin_starts(j) & ...
            aligned_spikes_cell{i} < bin_ends(j));
    end
end

% --- Final Outputs ---
aligned_spikes = vertcat(aligned_spikes_cell{:});

end