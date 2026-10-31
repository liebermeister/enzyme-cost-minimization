# Changelog

Version numbers have the form MAJOR.MINOR.PATCH. MAJOR changes when existing calls may stop working.

## 1.0.0 (not released yet)

First numbered version.

- Version function `ecm_version`
- `load_model_and_data_sbtab` and `save_model_and_data_sbtab` moved here from Parameter balancing for Matlab. The validation data file is now optional, and `write_delta_G0` is a real argument of `save_model_and_data_sbtab`
- `ecm_save_result_sbtab` accepts fluxes given as a vector
