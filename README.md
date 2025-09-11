# Analysis Pipeline for Tokens Task

## Overview

This repository contains the MATLAB analysis pipeline for the 'tokens' task. It is designed to process pre-consolidated neural and behavioral data (`session_data.mat`) to test hypotheses related to reward prediction error (RPE) and sensory prediction error (SPE) encoding in SC, SNc, and VTA.

## Getting Started

Before contributing to this repository, please review the developer guidelines in `code/AGENTS.md`. This document contains critical information regarding coding standards, data structure conventions, and project-specific practices.

## Repository Structure

-   **/code**: Contains all MATLAB source code.
    -   `run_tokens_analysis.m`: The main entry point for the analysis pipeline.
    -   **/utils**: General-purpose utility functions that are called by the main pipeline.
    -   **/examples**: Reference scripts that provide conceptual blueprints and stylistic patterns. These are not called directly by the pipeline.
-   **/config**: Contains configuration files.
    -   `session_manifest.csv`: Tracks the processing status of each recording session.
-   **/docs**: Contains project documentation.
    -   **/preprocessing_docs**: Documentation inherited from the upstream `neuro-preprocessing-pipeline`, including the critical data dictionary for `session_data.mat`.
-   **/figures**: Default output directory for generated figures. (Ignored by Git).

## Data

The analysis pipeline expects pre-consolidated data in the form of `session_data.mat` files. These files are the output of the separate `neuro-preprocessing-pipeline`. For a detailed description of the `session_data` structure and its fields, please consult the canonical data dictionary:

[**session_data Data Dictionary**](./docs/preprocessing_docs/session_data_dictionary.md)

## Pipeline Architecture

The pipeline is orchestrated by the main script `run_tokens_analysis.m`. It operates in several stages, with progress for each session tracked in the `config/session_manifest.csv` file. The script is idempotent, meaning it can be run multiple times without re-processing completed steps.

1.  **Neuron Screening:** For each session, selects task-modulated neurons using area-specific criteria.
    *   `screen_da_neurons.m`: Selects putative dopamine neurons based on firing rate and waveform shape.
    *   `screen_sc_neurons.m`: Selects task-modulated superior colliculus neurons and determines the recorded hemisphere.
2.  **Core Data Preparation:** The `prepare_core_data.m` script is called to generate analysis-ready data structures. It uses helper functions from the `/utils` directory to create event-aligned, binned spike rate matrices and processed pupil data for the selected neurons.
3.  **Define Task Conditions:** The `define_task_conditions.m` utility is used to create logical masks for various experimental conditions (e.g., trials with normal vs. uniform reward distributions).
4.  **Analysis Execution:** A suite of analyses are run on the prepared data:
    *   `analyze_baseline_comparison.m`: Compares baseline firing rates to post-event firing rates using ROC analysis.
    *   `analyze_roc_comparison.m`: Performs bin-by-bin ROC analysis to compare firing rates between key conditions (e.g., high vs. low RPE trials).
    *   `analyze_anova.m`: Performs a more complex N-way ANOVA to test for main effects and interactions of multiple task factors on firing rates.
5.  **Save Results:** The results of all analyses are saved to a single `analysis_results.mat` file in the `data/processed/{session_id}/` directory.

## Function Descriptions

This section provides a concise overview of the main `.m` files in the `/code` directory.

*   `run_tokens_analysis.m`: The main entry point for the analysis pipeline, which orchestrates neuron screening, data preparation, and analysis execution.
*   `screen_da_neurons.m`: Selects putative dopamine (DA) neurons based on firing rate and waveform characteristics.
*   `screen_sc_neurons.m`: Identifies task-modulated neurons in the superior colliculus (SC) and determines the recorded brain hemisphere.
*   `prepare_core_data.m`: A wrapper script that calls functions to prepare neuronal and pupil data for analysis.
*   `define_task_conditions.m`: Defines logical masks for various experimental conditions based on trial information and event times.
*   `analyze_baseline_comparison.m`: Compares baseline firing rates to post-event firing rates using ROC analysis.
*   `analyze_roc_comparison.m`: Performs a bin-by-bin ROC analysis to compare firing rates between two specified conditions.
*   `analyze_anova.m`: Performs an N-way ANOVA on firing rate data to test for main effects and interactions of task factors.
*   `generate_neuron_summary_pdf.m`: Generates a multi-page PDF with diagnostic plots for each individual neuron.
*   `plot_baseline_comparison.m`: Creates a summary figure visualizing the results of the baseline vs. post-event analysis.
*   `plot_roc_comparison.m`: Generates a summary figure visualizing the results of the between-condition ROC analysis.
*   `define_analysis_plan.m`: A configuration file that defines a plan for which analyses to run (currently unused).
*   `test_core_data_preparation.m`: A test script that verifies the functionality of the data preparation pipeline and generates diagnostic plots.
*   `test_neuron_diagnostics.m`: A test script that verifies the functionality of the neuron screening and diagnostic PDF generation workflow.

## Usage

1.  Ensure the `config/session_manifest.csv` is up to date.
2.  Run the main script from the MATLAB command window: `>> run_tokens_analysis`.
