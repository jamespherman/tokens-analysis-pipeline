# Aggregated Analysis Data Dictionary

This document defines the canonical data structure for the `aggregated_analysis_data.mat` file. This file is the output of the `aggregate_analysis_results.m` script and contains two top-level structs: `aggregated_sc_data` and `aggregated_snc_data`.

The structure within each of these is **dynamically generated**. The `aggregate_analysis_results.m` script recursively traverses the `.analysis` field of each individual `session_data.mat` file. The exact fields and analyses that are aggregated are determined by the `analysis_plan` defined in `define_task_conditions.m`. This makes the plan the single source of truth for what data appears in this aggregated file.

## Top-Level Structure (`aggregated_sc_data` and `aggregated_snc_data`)

| Field | Type | Description |
|---|---|---|
| `session_id` | `[N x 1] cell` | A cell array of strings indicating the `unique_id` of the session from which each neuron's data originated. [cite_start]`N` is the total number of neurons from all aggregated sessions for that brain area. [cite: 172-173] |
| `roc_comparison` | `struct` | Pooled results from the between-condition ROC analysis, as defined in `analysis_plan.roc_plan`. |
| `baseline_comparison` | `struct` | Pooled results from the baseline vs. post-event analysis, as defined in `analysis_plan.baseline_plan`. |
| `anova_results` | `struct` | Pooled results from the N-way ANOVA analysis, as defined in `analysis_plan.anova_plan`. |

---

### `roc_comparison` Structure

-   **Path:** `aggregated_sc_data.roc_comparison.(compName)`
-   **Description:** Contains the pooled results from `analyze_roc_comparison.m`. The specific comparisons run are defined in the `roc_plan` array within the `analysis_plan`. The `(compName)` in the path corresponds to the `.name` field from the plan (e.g., 'Dist_at_Cue').

| Field | Type | Description |
|---|---|---|
| `.sig` | `[N x T] double` | A matrix of significance values pooled from all sessions, where `N` is the total number of neurons and `T` is the number of time bins. |
| `.time_vector` | `[1 x T] double` | The canonical time vector for this analysis. [cite_start]The aggregation script uses the time vector from the first processed session that contains this analysis and stores it here, at the same level as the data matrix. [cite: 236-241] |
| `.cond_names` | `[1 x 2] cell` | The names of the two conditions that were compared (e.g., `{'is_normal_dist', 'is_uniform_dist'}`). This is taken from the first session containing this analysis. |

---

### `baseline_comparison` Structure

-   **Path:** `aggregated_sc_data.baseline_comparison.(eventName).(condName)`
-   **Description:** Contains the pooled results from `analyze_baseline_comparison.m`. The conditions analyzed are defined in `analysis_plan.baseline_plan`.

| Field | Type | Description |
|---|---|---|
| `.sig` | `[N x T] double` | A matrix of significance values pooled from all sessions. |
| `.time_vector` | `[1 x T] double`| The canonical time vector for this analysis, stored at the same level as the data matrix. [cite: 236-241] |

---

### `anova_results` Structure

-   **Path:** `aggregated_sc_data.anova_results.(eventName)`
-   **Description:** Contains the pooled results from `analyze_anova.m`. The specific events analyzed are defined in the `anova_plan` within the `analysis_plan`.

| Field | Type | Description |
|---|---|---|
| `.(p_value_name)` | `[N x T] double` | A matrix of p-values for a specific ANOVA term (e.g., `.p_rpe`). The specific p-value fields aggregated are defined in `analysis_plan.anova_plan.fields_to_aggregate`. |
| `.time_vector` | `[1 x T] double` | The canonical time vector for this analysis, stored at the same level as the p-value matrices. [cite: 236-241] |