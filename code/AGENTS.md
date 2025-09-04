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
