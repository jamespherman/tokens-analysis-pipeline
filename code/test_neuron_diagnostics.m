%% test_neuron_diagnostics.m
% This script serves as a step-by-step test of the neuron screening and
% diagnostic PDF workflow. It runs the necessary functions for a single,
% hardcoded session.

% Start timer
tic;

% In-line function to report timing
giveFeed = @(x)disp([num2str(toc) 's - ' x]);

%% Setup Paths
% Add the 'utils' directory to the path so that helper functions can be
% found.
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Setup
unique_id = 'Feynman_08_12_2025_SNc'; % Hardcoded session
giveFeed(sprintf('Testing diagnostic workflow for session: %s', ...
    unique_id));

%% Step 1: Load session data
giveFeed('Step 1: Loading session data...');
one_drive_path = findOneDrive;
session_data_path = fullfile(one_drive_path, 'Neuronal Data Analysis', ...
    unique_id, [unique_id '_session_data.mat']);

% Note: The following line will error if the data file does not exist.
% This is expected in a test environment without the actual data.
try
    data = load(session_data_path);
    session_data = data.session_data;
    giveFeed('Session data loaded successfully.');
catch ME
    warning('Could not load session data. Using dummy data instead.');
    warning(ME.message);
    % Create a dummy session_data struct for testing purposes
    session_data = struct('spikes', rand(100,1), 'neuron_ids', 1:5);
end

%% Step 2: Run neuron screening
giveFeed('Step 2: Running neuron screening...');
selected_neurons = []; % Initialize empty

if contains(unique_id, 'SNc')
    giveFeed('Session is SNc type. Running screen_da_neurons...');
    % Assumes screen_da_neurons is in the path
    selected_neurons = screen_da_neurons(session_data, unique_id);
    giveFeed('DA neuron screening complete.');
elseif contains(unique_id, 'SC')
    giveFeed(['Session is SC type. Loading gsac_data and running ' ...
        'screen_sc_neurons...']);
    % This part requires gsac_data, which is not available in this test.
    % We will simulate a placeholder for gsac_data.
    gsac_data_path = fullfile(one_drive_path, 'Neuronal Data Analysis', ...
        unique_id, [unique_id '_gsac_data.mat']);
    try
        gsac_data = load(gsac_data_path);
        giveFeed('gSac_jph task data loaded successfully.');
    catch ME_gsac
        warning('Could not load gsac_data. Using dummy data instead.');
        gsac_data = struct(); % Dummy gsac_data
    end

    % Assumes screen_sc_neurons is in the path
    selected_neurons = screen_sc_neurons(session_data, gsac_data);
    giveFeed('SC neuron screening complete.');
else
    error(['Unknown session type in unique_id. Cannot determine which ' ...
        'screening function to run.']);
end

%% Step 3: Calculate and Store Diagnostic Metrics
giveFeed('Step 3: Calculating and storing diagnostic metrics...');

% Calculate baseline firing rates and add to session_data
session_data.metrics.baseline_frs = calculate_baseline_fr(session_data);

% Calculate waveform metrics and add to session_data
nClusters = height(session_data.spikes.cluster_info);
for i_cluster = 1:nClusters
    % Note: Assuming a sampling rate of 30000 Hz, which is typical for this data
    session_data.metrics.wf_metrics(i_cluster) = ...
        calculate_waveform_metrics(session_data.spikes.wfMeans{i_cluster}, ...
        30000);
end
giveFeed('Diagnostic metrics calculated and stored.');


%% Step 4: Generate diagnostic PDF
giveFeed('Step 4: Generating diagnostic PDF...');
try
    generate_neuron_summary_pdf(session_data, selected_neurons, unique_id);
    giveFeed('Diagnostic PDF generation step complete.');
catch ME_pdf
    giveFeed('Error during PDF generation step.');
    rethrow(ME_pdf);
end

%% Done
giveFeed('Test complete. Check the generated PDF for session diagnostics.');
