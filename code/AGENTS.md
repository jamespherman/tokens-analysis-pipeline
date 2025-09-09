# Agent Instructions

## IMPORTANT: Data Structures

**Before writing or modifying any code that interacts with `session_data.mat` files, you MUST consult the data dictionary.**

The canonical definition for the main data structure, `session_data`, is located in:
`docs/preprocessing_docs/session_data_dictionary.md`

Do not make assumptions about the fields or layout of this structure. The documentation is the single source of truth. For example:
- The number of clusters/neurons should be derived from `session_data.spikes.cluster_info`, not a field like `nClusters`.
- Trial-related event times are in the `session_data.eventTimes` struct, not a `trials` struct.
- Spike times are stored as a single vector (`spikes.times`) and mapped to clusters via `spikes.clusters`, not as a cell array per neuron.

Consulting the documentation first will prevent bugs and ensure your code is compatible with the project's data standards.

---

When creating new MATLAB scripts in this directory, please ensure you add the `utils` directory to the MATLAB path. This is necessary for helper functions to be found.

You can use the following code snippet at the beginning of your script to achieve this:

```matlab
%% Setup Paths
% Add the 'utils' directory to the path so that helper functions can be
% found.
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
```
This will ensure that the path is added correctly, regardless of the current working directory.

---

## Coding Style

### Line Length
Each line of code should be limited to a maximum of 75 characters. For longer lines, please wrap them to a new line using `...`.

### Commenting
It is strongly preferred that comments are added on the line *above* the code to which they refer, rather than to the right of it.

---

## File Structure

All new MATLAB scripts should begin with a standard file header. This header should include a brief description of the script's purpose, the author's name, and the date of creation or last modification.

Example:
```matlab
%% my_script.m
%
% A brief one-line description of what the script does.
%
% A more detailed description can follow, explaining the inputs, outputs,
% and any key assumptions or dependencies.
%
% Author: Your Name
% Date: YYYY-MM-DD
%
```

Code should be organized into logical sections using `%%` headers. For example:

```matlab
%% Setup
% ... setup code here ...

%% Main Analysis
% ... main analysis code here ...

%% Plotting
% ... plotting code here ...
```

---
## Conventions

### **Analysis Pipeline**

#### **1. Selecting Trials by Task**
It is crucial to filter trials based on the specific task being analyzed. Neural responses are highly task-dependent, and calculations must be performed only on the relevant behavioral data. Use the `session_data.trialInfo.taskCode` field to select trials corresponding to a given task.

First, load the task code definitions, then create a logical mask to select the desired trial indices.

Example:
```matlab
% Load the structure containing all task codes
codes = initCodes();

% Find the indices of trials that belong to the 'tokens' task
tokens_trial_indices = session_data.trialInfo.taskCode == codes.uniqueTaskCode_tokens;

% Subsequent analyses should use these indices to filter the data
valid_cue_on_times = session_data.eventTimes.CUE_ON(tokens_trial_indices);
```

#### **2. Condition Mask Compatibility**
When creating logical masks for different experimental conditions, ensure they are compatible with the data they are intended to select. This means the masks should have the same number of elements as the trials for the specific task being analyzed.

For example, the `define_task_conditions.m` function is designed to generate condition masks specifically for 'tokens' trials. It first filters the session data to include only rewarded tokens trials and then creates masks that are the same length as this filtered set of trials. This ensures that the masks can be directly applied to `core_data` arrays, which are also filtered for tokens trials.

### **Interpreting Specific Task Data**
The `gSac_jph` task has special properties that can be leveraged during analysis. Memory-guided saccade trials within this task are intentionally placed at the neuron's estimated receptive/movement field center. This experimental design allows for two major simplifications:
1. The recorded side of the SC (`scSide`) can be inferred directly from the target location (e.g., targets in the left visual field imply a right SC recording).
2. All of these trials can be considered "contralateral" for the purposes of analysis, as they are designed to elicit a strong response from the recorded neurons.

---

## Project Root Directory

When writing code that needs to access files or directories relative to the project root, it is important to construct paths correctly. The project root for `tokens-analysis-pipeline` can be found by using the following MATLAB code:

```matlab
project_root = fullfile(findOneDrive, 'Code', 'tokens-analysis-pipeline');
```

This will ensure that paths are resolved correctly, regardless of the directory from which the code is executed.

For example, to access the `figures` directory from within a script, you should use:

```matlab
figures_dir = fullfile(project_root, 'figures');
```

This is especially important for scripts located in the `code/` directory, which are executed from `/code` in the runtime environment.

---

## Screening SC Neurons

When screening neurons from Superior Colliculus (SC) recordings, use the `screen_sc_neurons.m` function. This function has been designed to automate several key steps of the analysis for 'SC type' sessions.

### Function Usage
The function is called as follows:
```matlab
[selected_neurons, sig_epoch_comp, scSide] = screen_sc_neurons(session_data);
```

### Automatic `scSide` Determination
The function automatically determines the recorded side of the SC ('left' or 'right'). It does this by comparing the average visual-evoked firing rates for targets presented in the left versus the right visual field. The side with the stronger contralateral response is identified as the recorded hemisphere. The determined side is returned as the `scSide` output variable.

### Trial Selection
The function specifically uses memory-guided saccade trials from the `gSac_jph` task to perform neuron selection. It identifies these trials using the `taskCode` and the `isVisSac` field from `trialInfo`. If no such trials are found, it has a fallback mechanism to use any rewarded memory-guided saccade trials present in the session.

### `sig_epoch_comparison` Output
A key output of this function is `sig_epoch_comp`. This is a boolean matrix of size `[nClusters x 3]` that indicates whether a neuron showed a statistically significant change in firing rate for three critical comparisons:
1.  Visual epoch vs. Baseline
2.  Delay epoch vs. Baseline
3.  Saccade epoch vs. Baseline

This variable is crucial for the functional classification of SC neurons and should be stored for later analysis. The recommended practice is to store it in the `session_data` structure, for example:
```matlab
session_data.metadata.scSide = scSide;
session_data.metadata.sig_epoch_comparison = sig_epoch_comp;
```
This ensures that the results of the screening are saved with the session's data.

---

## Plotting Conventions

When creating complex, multi-panel figures for data verification or analysis, please adhere to the following conventions to ensure clarity, consistency, and scientific accuracy.

### 1. Never Smooth PSTHs
It is a sacred and holy and vitally important point that you **NEVER SMOOTH PSTHs**. Smoothing obscures event-locked timing information and defeats the purpose of using `barStairsFill`, which is designed to accurately represent the underlying binned spike counts.

- **Bad:** `movmean(psth, 5)`
- **Good:** `psth`

### 2. Use `barStairsFill` for PSTHs
For plotting Peristimulus Time Histograms (PSTHs), use the `barStairsFill.m` utility function. This function creates a staircase-style plot that accurately represents the binned nature of the data.

#### Plotting a Single PSTH
When plotting a single PSTH, the area under the curve should be filled to the zero line. This is the standard use case.

Example:
```matlab
% Plot a single PSTH with a fill to zero
barStairsFill(time_vector, psth, zeros(size(psth)));
```

#### Plotting Multiple PSTHs on the Same Axes
When plotting more than one PSTH on the same axes, **only the 'stairs' part should be visible, not the 'fill'**. This avoids obscuring the data. To achieve this, you must get the handles returned by the function and delete the fill and baseline objects.

Example:
```matlab
% Plotting two PSTHs (psth1 and psth2) on the same axes
hold on;

% Plot first PSTH as a line
h1 = barStairsFill(time_vector, zeros(size(psth1)), psth1);
delete(h1(1:2)); % Delete fill and baseline stairs
set(h1(3), 'Color', 'blue', 'LineWidth', 1.5); % Set color

% Plot second PSTH as a line
h2 = barStairsFill(time_vector, zeros(size(psth2)), psth2);
delete(h2(1:2)); % Delete fill and baseline stairs
set(h2(3), 'Color', 'red', 'LineWidth', 1.5); % Set color

legend({'Condition 1', 'Condition 2'});
```

### 3. De-clutter Axes in Multi-Panel Figures
To improve readability and reduce visual clutter in figures with multiple subplots arranged in a grid, remove redundant axis labels.
- For any given column of plots, only the bottom-most plot should have X-axis tick labels.
- For any given row of plots, only the left-most plot should have Y-axis tick labels.

Example:
```matlab
% h is an array of handles to the axes in a 2x3 grid
% Remove x-labels from the top row
set(h(1:3), 'XTickLabel', []);
% Remove y-labels from the middle and right columns
set(h([2,3,5,6]), 'YTickLabel', []);
```

### 4. Proportional Subplot Widths for Time-Series Data
When comparing different time epochs in adjacent subplots, the width of each subplot column should be proportional to the duration of the time window it represents. This provides an intuitive visual comparison of the temporal dynamics.

Do not use automated subplot tools like `subplot` or `mySubPlot` for this. Instead, manually calculate the `Position` vector `[left, bottom, width, height]` for each axis.

Example Logic:
```matlab
time_windows = [2, 2, 5.5]; % e.g., Cue (2s), Outcome (2s), Reward (5.5s)
proportions = time_windows / sum(time_windows);

% Calculate total available plotting width, accounting for margins/gaps
plot_area_width = 1 - left_margin - right_margin - (n_cols-1)*h_gap;

% Calculate individual column widths
col_widths = plot_area_width * proportions;

% Construct Position vectors using these widths
pos1 = [left_margin, bottom, col_widths(1), height];
h1 = axes('Position', pos1);
```

### 5. Use a Single, Informative `sgtitle`
Avoid using individual `title()` calls for each subplot. Instead, use a single, comprehensive main title for the entire figure using `sgtitle()`. This title should provide the key takeaway or context for the figure. Column and row labels can be added with `ylabel` on the leftmost plots and `xlabel` on the bottom plots.

### 6. Set Interpreter to 'none' for Titles with File Paths
When including filenames, unique IDs, or other text with special characters (like underscores `_`) in a title, always set the `'Interpreter'` property to `'none'`. This prevents MATLAB from treating the underscores as subscript indicators and ensures the text is displayed literally.

Example:
```matlab
unique_id = 'My_Session_2025_09_09';
sgtitle(sprintf('Data Verification for: %s', unique_id), 'Interpreter', 'none');
```
