# Analysis Pipeline for Tokens Task

## Overview

This repository contains the MATLAB analysis pipeline for the 'tokens' task. It is designed to process pre-consolidated neural and behavioral data (`session_data.mat`) to test hypotheses related to reward prediction error (RPE) and sensory prediction error (SPE) encoding in SC and SNc.

## Getting Started

Before contributing to this repository, please review the developer guidelines in `code/AGENTS.md`. This document contains critical information regarding coding standards, data structure conventions, and project-specific practices.

## Repository Structure

-   **/code**: Contains all MATLAB source code for the analysis.
    -   `run_tokens_analysis.m`: Main script for per-session processing.
    -   `aggregate_analysis_results.m`: Script to pool per-session results into a single dataset.
    -   `run_plotting_pipeline.m`: Main script for generating aggregated summary figures.
    -   **/utils**: General-purpose utility functions.
    -   **/examples**: Reference scripts that provide conceptual blueprints.
-   **/config**: Contains configuration files.
    -   `session_manifest.csv`: Tracks the processing status of each recording session.
-   **/data**: Contains processed data and, potentially, manually saved figures.
    -   **/processed**: Default output directory for aggregated data files.
    -   **/figures**: May contain manually saved summary plots.
-   **/docs**: Contains project documentation.
    -   **/preprocessing_docs**: Documentation for the upstream `neuro-preprocessing-pipeline`.
-   **/figures**: Default output directory for automatically generated per-session diagnostic figures. (Ignored by Git).

## Data and Figure Storage

This project generates several types of data and figures, which are stored in different locations.

-   **Per-Session Analysis Data**: Analysis results for each session are **not** stored within this repository. Instead, the main pipeline script (`run_tokens_analysis.m`) modifies the original `session_data.mat` files in-place, adding an `.analysis` field to the main struct. These files are located on a shared drive, typically found via the `findOneDrive()` utility function.

-   **Aggregated Analysis Data**: The results from all individual sessions are pooled into a single file located at:
    -   `data/processed/aggregated_analysis_data.mat`

-   **Per-Neuron Diagnostic Figures**: For each session, a PDF containing diagnostic plots for every neuron is automatically generated and saved to a session-specific subdirectory within the root `figures/` directory.
    -   Example: `figures/Feynman_08_05_2025_SNc/`

-   **Reward Distribution Figures**: For each session, a PDF visualizing the reward distributions used to define trial conditions is saved to the root `figures/` directory.
    -   Example: `figures/Feynman_08_05_2025_SNc_reward_distributions.pdf`

-   **Aggregated Summary Figures**: The pipeline can generate publication-quality summary plots (e.g., comparing SC vs. SNc). These figures are **not saved automatically**. The `run_plotting_pipeline.m` script will generate and display them, but the user must save them manually from the MATLAB figure window (e.g., to `data/figures/`).

## Pipeline Architecture

The analysis is a multi-stage process orchestrated by three main scripts. Progress for the first stage is tracked in `config/session_manifest.csv`.

1.  **Per-Session Processing (`run_tokens_analysis.m`)**:
    -   Iterates through each session in the manifest.
    -   **Neuron Screening**: Selects task-modulated neurons (`screen_da_neurons.m`, `screen_sc_neurons.m`).
    -   **Diagnostic PDF Generation**: Creates per-neuron summary plots (`generate_neuron_summary_pdf.m`).
    -   **Core Data Preparation**: Generates analysis-ready data structures (`prepare_core_data.m`).
    -   **Analysis Execution**: Runs a suite of analyses (ROC, ANOVA) on the prepared data.
    -   **Save Results**: Saves the analysis results by updating the session's original `session_data.mat` file in-place.

2.  **Results Aggregation (`aggregate_analysis_results.m`)**:
    -   Scans the manifest to find all sessions marked as 'complete'.
    -   Loads the `session_data.mat` file for each completed session.
    -   Extracts the `.analysis` field from each file.
    -   Pools the results by brain area (SC, SNc) into two large structs.
    -   Saves the aggregated data to `data/processed/aggregated_analysis_data.mat`.

3.  **Aggregated Plotting (`run_plotting_pipeline.m`)**:
    -   Loads the `aggregated_analysis_data.mat` file.
    -   Calls various `plot_aggregated_*.m` functions to generate summary figures.
    -   **Note**: Figures are displayed but not saved automatically.

## Function Descriptions

This section provides a concise overview of the main `.m` files in the `/code` directory.

*   `run_tokens_analysis.m`: The main entry point for the **per-session** analysis pipeline.
*   `aggregate_analysis_results.m`: Pools results from all sessions into a single aggregated data file.
*   `run_plotting_pipeline.m`: Generates and displays final summary figures from the aggregated data.
*   `screen_da_neurons.m`: Selects putative dopamine (DA) neurons.
*   `screen_sc_neurons.m`: Identifies task-modulated neurons in the superior colliculus (SC).
*   `prepare_core_data.m`: Prepares neuronal and pupil data for analysis.
*   `define_task_conditions.m`: Defines logical masks for experimental conditions.
*   `analyze_roc_comparison.m`: Performs bin-by-bin ROC analysis between conditions.
*   `analyze_anova.m`: Performs N-way ANOVA on firing rate data.
*   `generate_neuron_summary_pdf.m`: Generates a PDF of diagnostic plots for each neuron.
*   `plot_aggregated_roc_comparison.m`: Creates a summary figure from the aggregated ROC analysis results.
*   `plot_aggregated_anova.m`: Creates a summary figure from the aggregated ANOVA results.
*   `define_analysis_plan.m`: A configuration file that defines which analyses to run (currently unused).

## Usage

To run the full analysis from start to finish, execute the main scripts in the following order from the MATLAB command window:

1.  **Run Per-Session Analysis**:
    ```matlab
    >> run_tokens_analysis
    ```
    This will populate the `session_data.mat` files with analysis results. Monitor the command window for progress.

2.  **Aggregate Results**:
    ```matlab
    >> aggregate_analysis_results
    ```
    This creates the `aggregated_analysis_data.mat` file needed for plotting.

3.  **Generate Summary Plots**:
    ```matlab
    >> run_plotting_pipeline
    ```
    This will open the final summary figures. Remember to save them manually if needed.
