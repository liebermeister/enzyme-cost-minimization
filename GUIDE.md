# How do I ...?

A short guide to the most common tasks. Type `help` and a function name for details. The demo scripts are in the folder `enzyme-cost-minimization/demo`.

## First steps

```matlab
filenames = filenames_demo_ecm_simple;
options = struct('actions', 'ecm', 'ecm_scores', {{'emc4cm'}}, 'replace_cofactors', {{}}, 'fluctuations_safety_margin', 0);
[report, errors] = ecm_simple(filenames.model_data_file, filenames.validation_data_file, tempdir, options)
```

The input is one SBtab file with model and data. The results are written to SBtab files in the output folder.

## Run ECM

| Task | Function |
|---|---|
| Run parameter balancing or ECM on one SBtab file | `ecm_simple` |
| Run parameter balancing on a model in Matlab | `ecm_parameter_balancing` |
| Run ECM on a model in Matlab | `ecm_enzyme_cost_minimization` |

## Choose options and cost scores

| Task | Function |
|---|---|
| See all options and their default values | `ecm_default_options` |
| Complete the options | `ecm_update_options` |
| Read or save options | `ecm_options_read`, `ecm_options_save` |
| See the cost scores | `ecm_get_score` |

## Files and graphics

| Task | Function |
|---|---|
| Load a model and data from SBtab files | `load_model_and_data_sbtab` |
| Save a model and data to SBtab files | `save_model_and_data_sbtab` |
| Save results | `ecm_save_result_sbtab` |
| Read saved results | `ecm_load_result_sbtab` |
| Show graphics of results | `ecm_display` |
