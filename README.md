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

The pipeline is orchestrated by the main script `run_tokens_analysis.m`. It operates in several stages, with progress tracked in the `session_manifest.csv` file:

1.  **Neuron Screening:** Selects high-quality, task-modulated neurons based on area-specific criteria.
2.  **Data Preparation:** Generates event-aligned, binned spike rate matrices for the selected neurons.
3.  **Analysis Execution:** Runs a suite of analyses (e.g., discriminability, classification, CCA) on the prepared data.
4.  **Aggregation & Visualization:** Aggregates results across sessions and generates summary figures.

## Usage

1.  Ensure the `config/session_manifest.csv` is up to date.
2.  Run the main script from the MATLAB command window: `>> run_tokens_analysis`.
