# Changelog

Version numbers have the form MAJOR.MINOR.PATCH. MAJOR changes when existing calls may stop working.

## 1.0.0 (not released yet)

- `ecm_save_result_sbtab`: the model file contains the table of kinetic constants (`Parameter`) again (option `modular_rate_law_table` of mnt `network_to_sbtab`, whose default had changed to 0); this also affects the result files of model balancing (`balanced_model.tsv`), which lacked this table since the move of `network_to_sbtab` to mnt
- Standard concentration of thermodynamic quantities: new Config option `standard_concentration` (default `1 mM`, as before), used by `ecm_save_result_sbtab` (option `standard_concentration`) for the equilibrium constants and standard chemical potentials of the result files, declared by the table attribute `StandardConcentration` (also in the table of reaction Gibbs free energies and in the `keqs` table of `ecm_save_model_and_data_neos`). Needs mnt with the option `standard_concentration` of `network_to_sbtab`
- `load_model_and_data_sbtab`: data tables (`FluxData`, `ConcentrationConstraint`, `EnzymeCostWeight`, `MetaboliteConcentrationData`, `EnzymeConcentrationData`) are matched to the model elements by ID (columns `!Reaction`, `!Compound`); the order of rows no longer matters (before, the rows had to be in the order of the model elements). Missing enzyme cost weights are set to 1. Without a separate validation data file, validation data are read from the model document. Needs sbtab with `sbtab_table_get_values_by_id`
- Default options in the SBtab Config file `enzyme-cost-minimization/enzyme-cost-minimization_Config.tsv` (TableID `ConfigureECM`): `ecm_default_options` reads the default values from this file (on top of the defaults of parameter balancing) and only describes the options; options that depend on the network are still set by the function. Needs matlab-utils with `load_config_table` (later than 1.0.0); the dependency check warns if it is missing.
- New file `enzyme-cost-minimization/enzyme-cost-minimization_ToolInputOutput.tsv` (SBtab 1.2 table type ToolInputOutput) describing the input and output tables of ECM (for information only; not used by the code)
First numbered version.

- Example files with `e_coli` or `ecoli` in their names renamed to `escherichia-coli` (`escherichia_coli` in Matlab function and script names)
- Version function `ecm_version`
- `ecm_BASEDIR` and `ecm_RESOURCEDIR` renamed to `ecm_basedir` and `ecm_resourcedir`; they return the main directory and the `resources` directory, without a final file separator
- Folders `resources/model-files` and `resources/data-files` renamed to `resources/models` and `resources/data`
- Installation notes moved from `INSTALLATION` and the inner `README.md` to the main `README.md`; new `GUIDE.md`
- `CITATION.cff` and `.gitignore` added
- The main functions check the dependencies automatically, once per Matlab session
- New folder `tests` with tests (run with `runtests('tests')`)
- `load_model_and_data_sbtab` and `save_model_and_data_sbtab` moved here from Parameter balancing for Matlab. The validation data file is now optional, and `write_delta_G0` is a real argument of `save_model_and_data_sbtab`
- `ecm_save_result_sbtab` accepts fluxes given as a vector
