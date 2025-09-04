# Agent Instructions

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
