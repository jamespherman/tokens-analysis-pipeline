# Data Dictionary: `aggregated_analysis_data.mat`

This document provides a definitive guide to the structure of the `aggregated_analysis_data.mat` file. This file contains the pooled (aggregated) results from the analysis of individual recording sessions.

## Top-Level Structure

The `.mat` file contains two top-level variables:

1.  `aggregated_sc_data`: A struct containing all aggregated data for neurons identified as being in the Superior Colliculus (SC).
2.  `aggregated_snc_data`: A struct containing all aggregated data for neurons identified as being in the Substantia Nigra pars compacta (SNc).

Both of these structs share the exact same internal structure, as described below.

## Main Data Struct Structure (`aggregated_sc_data` and `aggregated_snc_data`)

Each of the top-level data structs is organized as follows. The fields are created dynamically based on the analyses that were run.

### Standard Fields

*   `session_id`: A `[N x 1]` cell array of strings, where `N` is the total number of neurons aggregated across all sessions. Each string is the unique identifier for the session from which the corresponding neuron's data was extracted. This allows tracing every data point back to its original session.

### Dynamic Analysis Fields

The other fields in the struct are named according to the specific analysis that was performed. For example, you might find fields such as:

*   `baseline_fr_hz`: An `[N x 1]` numeric array of baseline firing rates.
*   `roc_auc`: An `[N x C]` numeric array, where `C` is the number of conditions for a Receiver Operating Characteristic (ROC) analysis.
*   `psth`: A struct containing peristimulus time histogram (PSTH) data. The exact contents depend on the specific PSTH analysis run.

The key principle is that for any given field, the first dimension is always the neuron dimension, and it has `N` rows.

### The `time_vectors` Struct: A Special Case for Time-Series Data

This is a critical and often misunderstood part of the data structure.

Any analysis that produces a time-series output (like a PSTH or a windowed ROC analysis) will have its corresponding time vector stored in a special, nested struct called `time_vectors`.

*   **`time_vectors`**: This is a struct. The fields within `time_vectors` have the **same names** as the analysis fields they correspond to.

**Example:**

If you are working with a PSTH analysis stored in `aggregated_sc_data.psth_data`, the time vector for that analysis is **not** at the top level. Instead, it is located at:

```matlab
aggregated_sc_data.time_vectors.psth_data
```

This vector will have a size `[1 x T]`, where `T` is the number of time bins in the `psth_data` analysis `[N x T]`.

**Why this structure?**

Different analyses can have different time windows, binning, and alignment events. Storing time vectors this way ensures that every time-series data matrix is unambiguously associated with its correct time axis, preventing errors in plotting and interpretation.

**Golden Rule:** When plotting any data that has a time dimension, always look for its corresponding time vector in the `time_vectors` struct. Do not assume a single, global time vector exists.
