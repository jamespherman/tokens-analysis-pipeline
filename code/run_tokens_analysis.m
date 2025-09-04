%% run_tokens_analysis.m
%
% Main script to run the entire tokens task analysis pipeline.
%
% Author: Your Name
% Date: 2024-09-03
%

%% Setup
% Add necessary paths
addpath(genpath('code'));
addpath('config');

% Load the session manifest
manifest_path = 'config/session_manifest.csv';
manifest = readtable(manifest_path);

%% Iterate through sessions
for i = 1:height(manifest)

    % Check if the session needs processing
    if strcmp(manifest.screening_status{i}, 'pending')

        fprintf('Processing session: %s\n', manifest.session_id{i});

        % --- Placeholder: Load session_data.mat ---
        % session_data_path = fullfile(manifest.data_path{i}, manifest.session_id{i}, 'session_data.mat');
        % load(session_data_path);

        % --- Placeholder: Load gsac_data if needed ---
        % if contains(manifest.brain_area{i}, 'SC')
        %     gsac_data_path = fullfile(manifest.data_path{i}, manifest.session_id{i}, 'gsac_data.mat');
        %     load(gsac_data_path);
        % end

        % --- Neuron Screening ---
        if contains(manifest.brain_area{i}, 'SC')
            fprintf('--> Screening SC neurons...\n');
            % selected_neurons = screen_sc_neurons(session_data, gsac_data);
        elseif contains(manifest.brain_area{i}, 'SNc')
            fprintf('--> Screening DA neurons...\n');
            % selected_neurons = screen_da_neurons(session_data);
        else
            fprintf('--> No screening function defined for brain area: %s\n', manifest.brain_area{i});
            continue; % Skip to next session
        end

        % --- Placeholder: Save selected_neurons ---
        % output_dir = fullfile('data/interim', manifest.session_id{i});
        % if ~exist(output_dir, 'dir')
        %    mkdir(output_dir);
        % end
        % save(fullfile(output_dir, 'selected_neurons.mat'), 'selected_neurons');

        % --- Placeholder: Update manifest status ---
        % manifest.screening_status{i} = 'complete';
        % writetable(manifest, manifest_path);

        fprintf('--- Session %s processing complete ---\n\n', manifest.session_id{i});

    end
end

fprintf('All sessions processed.\n');
