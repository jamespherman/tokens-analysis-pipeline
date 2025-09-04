# Analysis Pipeline for Tokens Task

## Overview

This repository contains the MATLAB analysis pipeline for the 'tokens' task. It is designed to process pre-consolidated neural and behavioral data (`session_data.mat`) to test hypotheses related to reward prediction error (RPE) and sensory prediction error (SPE) encoding in SC, SNc, and VTA.

## Pipeline Architecture

The pipeline is orchestrated by the main script `run_tokens_analysis.m`. It operates in several stages, with progress tracked in the `session_manifest.csv` file:

1.  **Neuron Screening:** Selects high-quality, task-modulated neurons based on area-specific criteria.
2.  **Data Preparation:** Generates event-aligned, binned spike rate matrices for the selected neurons.
3.  **Analysis Execution:** Runs a suite of analyses (e.g., discriminability, classification, CCA) on the prepared data.
4.  **Aggregation & Visualization:** Aggregates results across sessions and generates summary figures.

## Usage

1.  Configure paths in `config/tokens_config.m`.
2.  Ensure the `config/session_manifest.csv` is up to date.
3.  Run the main script from the MATLAB command window: `>> run_tokens_analysis`.
